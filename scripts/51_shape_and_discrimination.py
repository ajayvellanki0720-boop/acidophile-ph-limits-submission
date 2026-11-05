"""Is the association linear across the pH range, and can it pick out an acidophile? -> DATA_DIR/shape_and_discrimination.csv
Every estimate so far is a single slope, which assumes the composite shifts the lower limit by the same amount everywhere from pH 0.4 to pH 9. Two reasons to doubt that. If the
mechanism is what the claim implies - fewer acidic surface residues where the surface meets the medium - the composite should matter most where the medium is actually acidic,
and the relationship should bend. And if it does not separate the extreme acidophiles from everything else, then whatever it measures is not about acid tolerance in the sense
the project set out to explain.
Three tests, all on the 32-protein exposed-periplasmic composite with genome-wide composition and data source controlled:
  1. SHAPE, by quintile of the composite rather than by band of the outcome. Binning on the outcome would restrict its range inside each bin and attenuate the slope by
     construction, so the predictor is binned instead and the mean limit is estimated per quintile by PGLS with the middle quintile as reference. This shows the shape without
     conditioning on what is being predicted.
  2. CURVATURE, as a squared term added to the linear model, which is the parametric version of the same question.
  3. DISCRIMINATION of extreme acidophiles (pH_min <= 3) by family-blocked cross-validated AUC, against the growth optimum and against genome-wide composition. AUC is used
     rather than R-squared because the question is a ranking one - given two species, does the composite put the more acid-tolerant one first - and because the acidophiles are
     a small minority, where R-squared is uninformative.
Usage: python3 scripts/51_shape_and_discrimination.py [n_repeats]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
z = lambda s: (s - s.mean()) / s.std()
cols = []
for r in MEM:
    if f"{r}|peri_exp|n" not in B or (r, "peri_exp") not in reflen:
        continue
    n = B[f"{r}|peri_exp|n"]
    cols.append(z(((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, "peri_exp")]))))
M = pd.concat(cols, axis=1)
B["x"] = M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
keep = (B.x.notna() & B.ph_min.notna()).to_numpy()
g = B[keep].reset_index(drop=True)
C = ph.vcv_fast(tree, list(B.ncbi_assembly))[np.ix_(keep, keep)]
y = g.ph_min.to_numpy(float)
SRCC = [f"src_{o}" for o in SRC]
use = [c for c in ["gw"] + SRCC if g[c].nunique() > 1]
base = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in use])
xz = ((g.x - g.x.mean()) / g.x.std()).to_numpy(float)
rows = []
print(f"{len(g)} bacteria with the composite on the tree; pH_min {y.min():.1f} to {y.max():.1f}\n")

print("1. mean lower limit by quintile of the composite (PGLS, middle quintile as reference)")
q = pd.qcut(g.x, 5, labels=False)
cen = 2
Xq = np.column_stack([base] + [(q == k).astype(float).to_numpy() for k in range(5) if k != cen])
fitq = ph.pgls(y, Xq, C)
names = [k for k in range(5) if k != cen]
lo = base.shape[1]
print(f"   quintile 3 (reference) n={int((q == cen).sum()):4d}  composite {g.x[q == cen].mean():+.3f}   mean pH_min {y[q == cen].mean():.2f}")
for j, k in enumerate(names):
    b, se, p = fitq["beta"][lo + j], fitq["se"][lo + j], fitq["p"][lo + j]
    rows.append({"test": "quintile", "group": f"Q{k + 1}", "n": int((q == k).sum()), "composite_mean": float(g.x[q == k].mean()),
                 "estimate": float(b), "se": float(se), "p": float(p)})
    print(f"   quintile {k + 1}             n={int((q == k).sum()):4d}  composite {g.x[q == k].mean():+.3f}   "
          f"{b:+.3f} +/- {se:.3f} pH vs reference (p {p:.4f}), raw mean pH_min {y[q == k].mean():.2f}")

print("\n2. curvature: a squared term added to the linear model")
fitl = ph.pgls(y, np.column_stack([base, xz]), C)
fitc = ph.pgls(y, np.column_stack([base, xz, xz ** 2]), C)
bl = fitl["beta"][-1]
b2, se2, p2 = fitc["beta"][-1], fitc["se"][-1], fitc["p"][-1]
rows.append({"test": "linear", "group": "composite", "n": len(g), "estimate": float(bl), "se": float(fitl["se"][-1]), "p": float(fitl["p"][-1])})
rows.append({"test": "quadratic", "group": "composite squared", "n": len(g), "estimate": float(b2), "se": float(se2), "p": float(p2)})
print(f"   linear only      {bl:+.4f} +/- {fitl['se'][-1]:.4f} per SD (p {fitl['p'][-1]:.2e})")
print(f"   with square      linear {fitc['beta'][-2]:+.4f}, square {b2:+.4f} +/- {se2:.4f} (p {p2:.3f})")

print("\n3. discrimination of extreme acidophiles (pH_min <= 3), family-blocked 10-fold cross-validated AUC")
d = g[g.family.notna() & g.ph_opt.notna()].reset_index(drop=True)
d["acid"] = (d.ph_min <= 3).astype(int)
print(f"   {len(d)} species with a family and an optimum, of which {int(d.acid.sum())} have pH_min <= 3 ({d.acid.mean():.1%})")
MODELS = {"composite alone": ["x"], "genome-wide composition alone": ["gw"], "growth optimum alone": ["ph_opt"],
          "optimum + composite": ["ph_opt", "x"], "optimum + composition": ["ph_opt", "gw"], "optimum + composite + composition": ["ph_opt", "x", "gw"]}


def auc(lab, sc):
    pos, neg = sc[lab == 1], sc[lab == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    r = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


fams = d.family.unique()
rng = np.random.default_rng(51)
res = {k: [] for k in MODELS}
for rep in range(NREP):
    assign = {f: i % 10 for i, f in enumerate(rng.permutation(fams))}
    fold = d.family.map(assign).to_numpy()
    for name, cs in MODELS.items():
        pred = np.full(len(d), np.nan)
        for k in range(10):
            tr, te = fold != k, fold == k
            if te.sum() < 5 or d.loc[tr, "acid"].nunique() < 2:
                continue
            Xtr, Xte = d.loc[tr, cs].to_numpy(float), d.loc[te, cs].to_numpy(float)
            mu, sg = Xtr.mean(0), Xtr.std(0); sg[sg == 0] = 1
            A = np.column_stack([np.ones(tr.sum()), (Xtr - mu) / sg])
            beta = np.linalg.lstsq(A, d.loc[tr, "acid"].to_numpy(float), rcond=None)[0]
            pred[te] = np.column_stack([np.ones(te.sum()), (Xte - mu) / sg]) @ beta
        # scores predict acidophily, which means LOW pH_min, so a lower limit is the positive class and the sign is already right
        m = ~np.isnan(pred)
        res[name].append(auc(d.acid.to_numpy()[m], pred[m]))
for name in MODELS:
    v = np.array(res[name])
    rows.append({"test": "AUC acidophile", "group": name, "n": len(d), "estimate": float(v.mean()), "se": float(v.std()), "p": np.nan})
    print(f"   {name:36s} AUC {v.mean():.3f} +/- {v.std():.3f}")
a, b_ = np.array(res["optimum + composite"]), np.array(res["growth optimum alone"])
print(f"   adding the composite to the optimum: {(a - b_).mean():+.4f} AUC, better in {int((a > b_).sum())}/{len(a)} fold assignments")
rows.append({"test": "AUC delta", "group": "composite over optimum", "n": len(d), "estimate": float((a - b_).mean()), "se": float((a - b_).std()),
             "p": np.nan, "improved_folds": int((a > b_).sum())})
pd.DataFrame(rows).to_csv(D / "shape_and_discrimination.csv", index=False)
print("\nwrote", D / "shape_and_discrimination.csv")
