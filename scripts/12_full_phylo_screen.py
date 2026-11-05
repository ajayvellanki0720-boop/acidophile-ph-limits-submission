"""Phylogeny-aware screen on the full marker trees: PGLS pH_min ~ pH_opt + source indicators (+ one feature), per domain, main tree,
two-sided p with Benjamini-Hochberg correction across features; Kdp tested one-sided (pre-specified direction) in Bacteria.
Usage: python3 scripts/12_full_phylo_screen.py [Bacteria] [Archaea]   (default: both; writes results/full/pgls_full_screen_<domain>.csv)"""
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
from acidophile.within_lineage import ALL_FEATURES

warnings.filterwarnings("ignore")
import os
D, T = ROOT / os.environ.get("DATA_DIR", "results/full"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_full")   # PHYLO_DIR=results/phylo_ml to use the ML tree
TAG = "" if T.name == "phylo_full" else "_" + T.name
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[df.ph_opt.notna()].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
SRC = ["src_thermobase", "src_acidb", "src_bacdive"]
FEATS = [f for f in ALL_FEATURES + ["gc"] if f in df.columns]
rows, kdp_rows = [], []
DOMS = sys.argv[1:] or ["Bacteria", "Archaea"]
for dom in DOMS:
    tree = Phylo.read(T / f"{dom}.nwk", "newick")
    tips = {t.name for t in tree.get_terminals()}
    g = df[(df.domain == dom) & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
    C = ph.vcv_from_tree(tree, list(g.ncbi_assembly))
    y = g.ph_min.to_numpy(float)
    covs = ["ph_opt"] + [s for s in SRC if g[s].nunique() > 1]
    base = ph.pgls(y, np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in covs]), C)
    print(f"{dom}: n={len(g)} | sources {g.origin.value_counts().to_dict()} | baseline lambda={base['lambda']:.2f} | beta(pH_opt)={base['beta'][1]:.2f}")
    for f in [f for f in FEATS if g[f].nunique() > 2]:
        z = ((g[f] - g[f].mean()) / g[f].std()).to_numpy(float)
        X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in covs] + [z])
        r = ph.pgls(y, X, C); k = X.shape[1] - 1
        rows.append({"domain": dom, "feature": f, "n": len(g), "beta_per_sd": r["beta"][k], "se": r["se"][k], "p_two_sided": r["p"][k],
                     "lambda": r["lambda"], "prevalence_gt0": float((g[f] > 0).mean())})
R = pd.DataFrame(rows)
for dom, idx in R.groupby("domain").groups.items():
    p = R.loc[idx, "p_two_sided"].to_numpy(); o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    R.loc[idx, "q_BH"] = np.minimum(q, 1)
for dom in DOMS:
    R[R.domain == dom].sort_values("p_two_sided").to_csv(D / f"pgls_full_screen_{dom}{TAG}.csv", index=False)
pd.set_option("display.width", 220)
for dom in DOMS:
    print(f"\n{dom}: features tested {int((R.domain == dom).sum())}"); print(R[R.domain == dom].sort_values("p_two_sided").head(8).round(3).to_string(index=False))
k = R[(R.domain == "Bacteria") & (R.feature == "hm_rate_kdp")]
if len(k):
    r = k.iloc[0]
    print(f"\nPre-specified Kdp (Bacteria, one-sided negative): beta={r.beta_per_sd:.3f} per SD, p={stats.t.cdf(r.beta_per_sd / r.se, r.n - 6):.4f}, prevalence {r.prevalence_gt0:.2f}")
