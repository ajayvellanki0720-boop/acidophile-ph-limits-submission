"""Diagnostics for the full-tree PGLS signal (Kdp negative, Trk/Ktr positive) in Bacteria.
(a) exact statistics, (b) leave-one-phylum-out stability, (c) permutation within families / phyla, (d) calibration against traits simulated
by Brownian motion on the same tree and mapped to each feature's own marginal distribution, (e) family-blocked prediction gain.
Usage: python3 scripts/13_kdp_trk_diagnostics.py [n_sims]"""
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
from acidophile.models import compare_models

warnings.filterwarnings("ignore")
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
import os
D, T = ROOT / os.environ.get("DATA_DIR", "results/full"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_full")   # PHYLO_DIR=results/phylo_ml to use the ML tree
TAG = "" if T.name == "phylo_full" else "_" + T.name
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[df.ph_opt.notna() & (df.domain == "Bacteria")].reset_index(drop=True)
tree = Phylo.read(T / "Bacteria.nwk", "newick")
df = df[df.ncbi_assembly.isin({t.name for t in tree.get_terminals()})].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
covs = ["ph_opt", "src_thermobase", "src_acidb", "src_bacdive"]
C = ph.vcv_from_tree(tree, list(df.ncbi_assembly))
y = df.ph_min.to_numpy(float)
rng = np.random.default_rng(11)


def zs(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std()


def design(z, d=df):
    return np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in covs] + [z])


FEATS = {"Kdp": "hm_rate_kdp", "Trk/Ktr": "hm_rate_trk_ktr"}
out = []
for label, f in FEATS.items():
    z = zs(df[f]); r = ph.pgls(y, design(z), C); k = design(z).shape[1] - 1
    out.append({"feature": label, "test": "PGLS main tree (n=%d)" % len(df), "beta_per_sd": round(r["beta"][k], 3), "p": f"{r['p'][k]:.2e}", "note": f"lambda={r['lambda']:.2f}"})
    # (b) leave one phylum out
    bs = []
    for phy in df.phylum.dropna().unique():
        keep = (df.phylum != phy).to_numpy()
        if keep.sum() < 40:
            continue
        sub = df[keep].reset_index(drop=True); Cs = C[np.ix_(keep, keep)]
        rr = ph.pgls(sub.ph_min.to_numpy(float), design(zs(sub[f]), sub), Cs); bs.append((phy, rr["beta"][-1], rr["p"][-1]))
    b = pd.DataFrame(bs, columns=["phylum", "beta", "p"])
    worst = b.loc[b.p.idxmax()]
    out.append({"feature": label, "test": "leave-one-phylum-out (%d fits)" % len(b), "beta_per_sd": f"{b.beta.min():.3f} to {b.beta.max():.3f}",
                "p": f"max p {b.p.max():.1e}", "note": f"weakest when dropping {worst.phylum}"})
    # (c) permutation of the feature within strata, using OLS partial association (residual of ph_min on covariates)
    Xc = np.column_stack([np.ones(len(df))] + [df[c].to_numpy(float) for c in covs])
    res = y - Xc @ np.linalg.lstsq(Xc, y, rcond=None)[0]
    x = df[f].to_numpy(float)
    obs = stats.spearmanr(res, x)[0]
    for strat in ("family", "phylum"):
        s = df[strat].fillna("NA").to_numpy(); idx = [np.where(s == k_)[0] for k_ in np.unique(s)]
        informative = sum(len(ix) > 1 and len(set(x[ix])) > 1 for ix in idx)
        null = np.empty(3000)
        for i in range(3000):
            xp = x.copy()
            for ix in idx:
                xp[ix] = rng.permutation(x[ix])
            null[i] = stats.spearmanr(res, xp)[0]
        out.append({"feature": label, "test": f"permute within {strat}", "beta_per_sd": f"rho {obs:.3f}", "p": f"{(np.abs(null) >= abs(obs)).mean():.4f}",
                    "note": f"{informative} {strat} strata have variation in the feature"})
    # (d) phylogenetic-simulation calibration
    L = np.linalg.cholesky(C / np.mean(np.diag(C)) + 1e-8 * np.eye(len(df)))
    srt = np.sort(x); ts = []
    for _ in range(NS):
        sim = L @ rng.normal(size=len(df))
        xs = srt[np.argsort(np.argsort(sim))]            # same marginal distribution as the real feature, tree-structured
        rr = ph.pgls(y, design(zs(xs)), C); ts.append(rr["t"][-1])
    ts = np.array(ts); tobs = r["t"][k]
    out.append({"feature": label, "test": f"calibration vs {NS} tree-simulated null features", "beta_per_sd": f"t_obs {tobs:.2f}",
                "p": f"{((np.abs(ts) >= abs(tobs)).sum() + 1) / (NS + 1):.4f}", "note": f"null |t| 95th pct {np.percentile(np.abs(ts), 95):.2f}; share of null |t|>1.96: {(np.abs(ts) > 1.96).mean():.2f}"})
# (e) prediction gain, family-blocked
import acidophile.models as m
base_cols = df.copy()
res0, _ = compare_models(base_cols, group_col="family", scheme="family_blocked")
out.append({"feature": "(both)", "test": "family-blocked CV: M1 (opt only) RMSE", "beta_per_sd": round(res0[res0.model == "M1_opt_only"].rmse.iloc[0], 3), "p": "", "note": f"n={len(df)}"})
OUTDF = pd.DataFrame(out)
OUTDF.to_csv(D / f"kdp_trk_diagnostics{TAG}.csv", index=False)
pd.set_option("display.width", 240); pd.set_option("display.max_colwidth", 90)
print(OUTDF.to_string(index=False))
