"""Ortholog-matched test: does the acidic-residue association hold WITHIN conserved protein families (same proteins in every genome), and is it stronger in membrane
than soluble families? PGLS (source indicators + family's acidic excess), per family, bacteria; plus random-effects pooling by class of family and the optimum-conditioned version.
Usage: python3 scripts/33b_ortholog_test.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
NAMES = {"PF00344": "SecY", "PF00119": "ATP synthase a", "PF00137": "ATP synthase c", "PF00115": "cytochrome oxidase I", "PF00361": "NADH-quinone TM subunits", "PF02096": "YidC",
         "PF00009": "EF-Tu/EF-G", "PF00154": "RecA", "PF00118": "chaperonin GroEL", "PF00006": "ATP synthase alpha/beta", "PF00012": "Hsp70/DnaK"}
MEM = ["PF00344", "PF00119", "PF00137", "PF00115", "PF00361", "PF02096"]; SOL = ["PF00009", "PF00154", "PF00118", "PF00006", "PF00012"]
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly").merge(pd.read_csv(D / "ortholog_charge.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float)
rows = []
for fam in MEM + SOL:
    col = f"{fam}_acidic_excess"
    if col not in B:
        continue
    for lab, with_opt in (("total effect", False), ("optimum as covariate", True)):
        m = B[col].notna().to_numpy() & (B.ph_opt.notna().to_numpy() if with_opt else True)
        if m.sum() < 200:
            continue
        g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
        cov = (["ph_opt"] if with_opt else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
        X = np.column_stack([np.ones(len(g))] + [g[c] for c in cov])
        b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, g[[col]].to_numpy(float))
        rows.append({"family": NAMES[fam], "pfam": fam, "kind": "membrane" if fam in MEM else "soluble", "model": lab, "n": len(g), "sd_pts": 100 * g[col].std(), "beta_per_sd": b[0], "se": se[0], "p": p[0]})
R = pd.DataFrame(rows); R.to_csv(D / "ortholog_test.csv", index=False)
pd.set_option("display.width", 200)
for lab in ("total effect", "optimum as covariate"):
    x = R[R.model == lab]; print(f"\n{lab}: acidic excess of the conserved family vs pH_min (beta per SD of that family)")
    print(x[["kind", "family", "n", "sd_pts", "beta_per_sd", "se", "p"]].round(4).to_string(index=False))
    for kind in ("membrane", "soluble"):
        y = x[x.kind == kind]; w = 1 / y.se ** 2; fe = (w * y.beta_per_sd).sum() / w.sum(); Q = (w * (y.beta_per_sd - fe) ** 2).sum(); k = len(y) - 1
        tau2 = max(0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (y.se ** 2 + tau2); mu = (wr * y.beta_per_sd).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
        print(f"   {kind:9s} families pooled: {mu:+.3f} ± {sem:.3f} per SD (p {2 * stats.norm.sf(abs(mu / sem)):.1e}); positive in {(y.beta_per_sd > 0).sum()}/{len(y)}; I² {max(0, (Q - k) / Q):.0%}")
