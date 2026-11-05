"""Panel exposure test (data from 36a) -> DATA_DIR/structure_panel_test.csv.
SCALE WARNING, correction 15. The predictor is built in percentage points but `beta_per_pt` holds a coefficient per STANDARD DEVIATION of that share, because
ph.pgls_scan standardises what it tests and nothing here divides the SD back out. The stored `sd_pts` column is the divisor: per point = beta_per_pt / sd_pts.
Per SD is the right scale for comparing classes and references, which is what this script does, so the values are correct and the column name is not.
PGLS of pH_min on the acidic excess (D+E-K-R per residue, percentage points) of each structural class of
each reference, in bacteria on the marker tree, with source indicators; three models: total effect, optimum as covariate, and total effect plus genome-wide acidic excess.
Then random-effects pooling (DerSimonian-Laird) of each class across the seven membrane references, Benjamini-Hochberg over the classes tested, and the Escherichia vs Bacillus
SecY comparison as a reference-species control. A class is scored in a genome only if at least max(8, half the reference positions of that class) align.
Usage: python3 scripts/36b_structure_panel_test.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
SOL = ["EF-Tu", "GroEL", "RecA"]
MCL = ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"]
SCL = ["surface", "buried"]
reflen = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
rows = []
for ref, classes in [(r, MCL) for r in MEM] + [(r, SCL) for r in SOL]:
    for c in classes:
        if f"{ref}|{c}|n" not in B or (ref, c) not in reflen:
            continue
        n = B[f"{ref}|{c}|n"]
        need = max(8, 0.5 * reflen[(ref, c)])
        x = ((B[f"{ref}|{c}|D"] + B[f"{ref}|{c}|E"] - B[f"{ref}|{c}|K"] - B[f"{ref}|{c}|R"]) / n).where(n >= need) * 100
        for lab in ("total effect", "optimum as covariate", "plus genome-wide acidic excess"):
            m = x.notna().to_numpy() & (B.ph_opt.notna().to_numpy() if lab == "optimum as covariate" else True)
            if m.sum() < 200:
                continue
            g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
            cov = (["ph_opt"] if lab == "optimum as covariate" else []) + (["gw"] if lab.startswith("plus") else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
            X = np.column_stack([np.ones(len(g))] + [g[k] for k in cov])
            xv = x[m].to_numpy(float)
            b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, xv[:, None])
            rows.append({"ref": ref, "kind": "membrane" if ref in MEM else "soluble", "class": c, "model": lab, "n": len(g),
                         "ref_positions": int(reflen[(ref, c)]), "mean_aligned": float(n[m].mean()), "sd_pts": xv.std(), "beta_per_pt": b[0], "se": se[0], "p": p[0]})
R = pd.DataFrame(rows)
for lab, gr in R.groupby("model"):
    R.loc[gr.index, "q"] = stats.false_discovery_control(gr.p.to_numpy(), method="bh")
R.to_csv(D / "structure_panel_test.csv", index=False)
pd.set_option("display.width", 220)


def pool(y):
    w = 1 / y.se ** 2; fe = (w * y.beta_per_pt).sum() / w.sum(); Q = (w * (y.beta_per_pt - fe) ** 2).sum(); k = len(y) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (y.se ** 2 + tau2); mu = (wr * y.beta_per_pt).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), (max(0.0, (Q - k) / Q) if k > 0 and Q > 0 else 0.0)


for lab in ("total effect", "optimum as covariate", "plus genome-wide acidic excess"):
    x = R[R.model == lab]
    print(f"\n=== {lab}: pH_min vs acidic excess per percentage point")
    print(x[["kind", "ref", "class", "n", "ref_positions", "mean_aligned", "sd_pts", "beta_per_pt", "se", "p", "q"]].round(4).to_string(index=False))
    print("  pooled across the membrane references (random effects):")
    for c in MCL:
        y = x[(x.kind == "membrane") & (x["class"] == c)]
        if len(y) < 2:
            continue
        mu, sem, pv, i2 = pool(y)
        print(f"    {c:9s} {mu:+.3f} +/- {sem:.3f} (p {pv:.1e}); positive in {(y.beta_per_pt > 0).sum()}/{len(y)}; I2 {i2:.0%}")
    y = x[x.kind == "soluble"]
    for c in SCL:
        z = y[y["class"] == c]
        if len(z) >= 2:
            mu, sem, pv, i2 = pool(z)
            print(f"    soluble {c:7s} {mu:+.3f} +/- {sem:.3f} (p {pv:.1e}); positive in {(z.beta_per_pt > 0).sum()}/{len(z)}; I2 {i2:.0%}")
    e, b = x[(x.ref == "SecY") & (x["class"] == "peri_exp")], x[(x.ref == "SecY_Bsub") & (x["class"] == "peri_exp")]
    if len(e) and len(b):
        d = e.beta_per_pt.iloc[0] - b.beta_per_pt.iloc[0]; sd = np.hypot(e.se.iloc[0], b.se.iloc[0])
        print(f"  SecY reference species, exposed periplasmic: Escherichia {e.beta_per_pt.iloc[0]:+.3f}+/-{e.se.iloc[0]:.3f} ({e.ref_positions.iloc[0]} positions) vs "
              f"Bacillus {b.beta_per_pt.iloc[0]:+.3f}+/-{b.se.iloc[0]:.3f} ({b.ref_positions.iloc[0]} positions); difference {d:+.3f}+/-{sd:.3f} (p {2 * stats.norm.sf(abs(d / sd)):.2f})")
