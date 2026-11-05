"""Robustness of the periplasmic-exposed signal (see 34a/34b): (1) joint PGLS of pH_min on all classes of one protein at once, in the genomes that score every class;
(2) the periplasmic-exposed association restricted to genomes aligning >= 80 % of that class, and by phylum group (Pseudomonadota alone; everything else); (3) optimum-conditioned
joint model. Usage: python3 scripts/34c_structure_exposure_robust.py"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
REFLEN = {("YidC", "peri_exp"): 158, ("SecY", "peri_exp"): 71}
df = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "structure_exposure.csv"), on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))


def ex(ref, c):
    n = B[f"{ref}|{c}|n"]
    return ((B[f"{ref}|{c}|D"] + B[f"{ref}|{c}|E"] - B[f"{ref}|{c}|K"] - B[f"{ref}|{c}|R"]) / n).where(n >= 15) * 100, n


def fit(m, cols, with_opt=False):
    m = m & (B.ph_opt.notna().to_numpy() if with_opt else True)
    g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
    cov = (["ph_opt"] if with_opt else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in cov])
    # joint model: first column of the scan is tested with the others as fixed covariates -> run each in turn
    out = {}
    for k, c in enumerate(cols):
        others = [x for j, x in enumerate(cols) if j != k]
        Xk = np.column_stack([X] + [g[o] for o in others]) if others else X
        b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), Xk, g[[c]].to_numpy(float))
        out[c] = (b[0], se[0], p[0])
    return len(g), out


print("(1) joint models, all classes of one protein entered together (beta per percentage point)")
for ref, classes in (("YidC", ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"]), ("SecY", ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"])):
    cols = []
    for c in classes:
        B[f"{ref}_{c}"], _ = ex(ref, c); cols.append(f"{ref}_{c}")
    m = B[cols].notna().all(axis=1).to_numpy()
    for lab, wo in (("total", False), ("optimum-conditioned", True)):
        n, out = fit(m, cols, wo)
        print(f"  {ref} {lab} (n={n}): " + " | ".join(f"{c.split('_',1)[1]} {v[0]:+.3f}±{v[1]:.3f} (p {v[2]:.1e})" for c, v in out.items()))

print("\n(2) periplasmic-exposed class, restricted samples (beta per percentage point, total effect)")
for ref in ("YidC", "SecY"):
    c = f"{ref}_peri_exp"; B[c], n = ex(ref, "peri_exp"); full = REFLEN[(ref, "peri_exp")]
    for lab, m in (("all scored genomes", B[c].notna()), ("aligned >=80% of the class", B[c].notna() & (n >= 0.8 * full)),
                   ("Pseudomonadota only", B[c].notna() & (B.phylum == "Pseudomonadota")), ("all other phyla", B[c].notna() & (B.phylum != "Pseudomonadota"))):
        m = m.to_numpy()
        if m.sum() < 150:
            print(f"  {ref} {lab}: n={m.sum()} (too few)"); continue
        n_, out = fit(m, [c]); b, se, p = out[c]
        print(f"  {ref} {lab}: n={n_}  {b:+.3f} ± {se:.3f} (p {p:.1e})")
