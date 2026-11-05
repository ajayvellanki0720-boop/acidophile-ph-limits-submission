"""Robustness of acidic excess / mean pI / acidic-residue fraction and Kdp in the broadened bacterial set: trimming the most acidic proteomes (halophile-like), dropping
halophile-named genera, and permutation within phyla / families. Usage: python3 scripts/31_acidic_robustness.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
C = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float)
FE = ["acidic_excess", "mean_pI", "frac_acidic", "hm_rate_kdp"]
rows = []


def run(label, mask):
    idx = np.where(np.asarray(mask))[0]; g = B.iloc[idx]
    cols = [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
    X = np.column_stack([np.ones(len(idx))] + [g[c] for c in cols]); b, se, p, lam = ph.pgls_scan(C[np.ix_(idx, idx)], Y[idx], X, g[FE].to_numpy(float))
    for f, bb, ss, pp in zip(FE, b, se, p):
        rows.append({"analysis": label, "n": len(idx), "feature": f, "beta_per_sd": bb, "se": ss, "p": pp})
    print(f"{label:70s} n={len(idx):4d} | " + " | ".join(f"{f.replace('hm_rate_', '')} {bb:+.3f} (p {pp:.0e})" for f, bb, pp in zip(FE, b, p)))


q = B.acidic_excess.quantile
run("all", np.ones(len(B), bool))
for cut in (0.99, 0.95, 0.90):
    run(f"drop most acidic proteomes (top {int(round((1 - cut) * 100))}% acidic excess)", B.acidic_excess <= q(cut))
run("drop top 10% and bottom 10% acidic excess", (B.acidic_excess <= q(.9)) & (B.acidic_excess >= q(.1)))
hal = B.genus.fillna("").str.contains("Salinibacter|Halo|Salini|Salicola|Chromohalobacter|Halomonas|Halobacillus|Natranaerobius|Salinivibrio|Halorhabdus", case=False)
run(f"drop halophile-named genera ({int(hal.sum())} species)", ~hal)
Xc = np.column_stack([np.ones(len(B))] + [B[f"src_{o}"] for o in SRC if B[f"src_{o}"].nunique() > 1]); res = Y - Xc @ np.linalg.lstsq(Xc, Y, rcond=None)[0]
rng = np.random.default_rng(3)
for f in ("acidic_excess", "mean_pI", "hm_rate_kdp"):
    for strat in ("phylum", "family"):
        s = B[strat].fillna("NA").to_numpy(); idx = [np.where(s == k)[0] for k in np.unique(s)]; x = B[f].to_numpy(float); obs = stats.spearmanr(res, x)[0]; null = np.empty(1500)
        for i in range(1500):
            xp = x.copy()
            for ix in idx:
                xp[ix] = rng.permutation(x[ix])
            null[i] = stats.spearmanr(res, xp)[0]
        pv = (np.abs(null) >= abs(obs)).mean()
        rows.append({"analysis": f"permute within {strat}", "n": len(B), "feature": f, "beta_per_sd": obs, "se": np.nan, "p": pv})
        print(f"   permute within {strat:6s} | {f:14s} rho {obs:+.3f} p {pv:.4f}{' (no permutation reached the observed value; p < 1/1500)' if pv == 0 else ''}")
pd.DataFrame(rows).to_csv(D / "acidic_robustness.csv", index=False)
