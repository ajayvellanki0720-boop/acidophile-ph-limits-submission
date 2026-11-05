"""Replication of the potassium-system leads (Kdp negative, Trk/Ktr positive) on the expanded set, with the earlier hypotheses fixed in advance.
 A. Joint PGLS (pH_opt, source indicators if >1 source, + Kdp, Trk/Ktr, His fraction) on: new species only (parser-extracted descriptions), the previously analysed
    literature species, database rows, and everything; also with genome size and GC added.
 B. Recurrence: the same model within each class / phylum with enough species, random-effects pooling, heterogeneity, sign counts.
 C. Calibration against tree-simulated null features on the new species (n_sims), and permutation within phyla / families.
Usage: python3 scripts/25_replication_broad.py [n_sims]   (env DATA_DIR=results/broad PHYLO_DIR=results/phylo_broad)"""
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
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 300
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
B["log_genome"] = np.log10(B.genome_size_bp)
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    B[f"src_{o}"] = (B.origin == o).astype(float)
z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)
FEATS = {"Kdp": "hm_rate_kdp", "Trk/Ktr": "hm_rate_trk_ktr", "His": "frac_His"}
rng = np.random.default_rng(31)


def fit(mask, extra=(), feats=FEATS, tree_=tree):
    g = B[mask].reset_index(drop=True)
    cov = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive", "src_literature_auto") if g[c].nunique() > 1] + list(extra)
    keep = [k for k, c in feats.items() if g[c].std() > 0]
    C = ph.vcv_fast(tree_, list(g.ncbi_assembly))
    X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in cov] + [z(g[feats[k]]) for k in keep])
    r = ph.pgls(g.ph_min.to_numpy(float), X, C); k0 = X.shape[1] - len(keep)
    return {k: (r["beta"][k0 + i], r["se"][k0 + i], r["p"][k0 + i]) for i, k in enumerate(keep)} | {"_n": len(g), "_lambda": r["lambda"]}, g, C


rows = []
subsets = {"new species only (parser-extracted)": B.origin == "literature_auto", "previously analysed literature species": B.origin == "literature",
           "database rows (ThermoBase, AciDB, BacDive)": B.origin.isin(["thermobase", "acidb", "bacdive"]), "all bacteria": np.ones(len(B), bool)}
for label, mask in subsets.items():
    for sens, extra in (("", ()), (" + genome size, GC", ("log_genome", "gc"))):
        res, g, C = fit(np.asarray(mask), extra)
        for k in ("Kdp", "Trk/Ktr", "His"):
            if k in res:
                b, se, p = res[k]; rows.append({"subset": label + sens, "n": res["_n"], "feature": k, "beta_per_sd": b, "se": se, "p_two_sided": p,
                                                "p_one_sided_expected_direction": stats.norm.cdf(b / se) if k == "Kdp" else (stats.norm.sf(b / se) if k == "Trk/Ktr" else np.nan), "lambda": res["_lambda"]})
A = pd.DataFrame(rows); A.to_csv(D / "replication_by_subset.csv", index=False)
pd.set_option("display.width", 230)
print("A. joint PGLS (Kdp negative / Trk-Ktr positive pre-specified in direction; His exploratory)"); print(A.round(5).to_string(index=False))

# C1. calibration on the new species
mask = (B.origin == "literature_auto").to_numpy(); res, g, C = fit(mask); y = g.ph_min.to_numpy(float)
cov = ["ph_opt"]; Xb = np.column_stack([np.ones(len(g)), g.ph_opt])
L = np.linalg.cholesky(C / np.mean(np.diag(C)) + 1e-8 * np.eye(len(g)))
cal = []
for k, c in FEATS.items():
    x = g[c].to_numpy(float)
    if x.std() == 0:
        continue
    srt = np.sort(x); ts = [ph.pgls(y, np.column_stack([Xb, z(srt[np.argsort(np.argsort(L @ rng.normal(size=len(g))))])]), C)["t"][-1] for _ in range(NS)]
    tobs = ph.pgls(y, np.column_stack([Xb, z(x)]), C)["t"][-1]
    cal.append({"feature": k, "t_obs": tobs, "empirical_p": ((np.abs(ts) >= abs(tobs)).sum() + 1) / (NS + 1), "null_95th_abs_t": np.percentile(np.abs(ts), 95), "n_sims": NS})
Cal = pd.DataFrame(cal); Cal.to_csv(D / "replication_calibration_new_species.csv", index=False)
print("\nC. calibration vs tree-simulated features (new species, single-feature model with pH_opt):"); print(Cal.round(4).to_string(index=False))
# C2. permutation within strata (OLS partial association on the new species)
Xc = np.column_stack([np.ones(len(g)), g.ph_opt]); resid = y - Xc @ np.linalg.lstsq(Xc, y, rcond=None)[0]
for strat in ("phylum", "family"):
    s = g[strat].fillna("NA").to_numpy(); idx = [np.where(s == k)[0] for k in np.unique(s)]
    for k, c in FEATS.items():
        x = g[c].to_numpy(float); obs = stats.spearmanr(resid, x)[0]; null = np.empty(2000)
        for i in range(2000):
            xp = x.copy()
            for ix in idx:
                xp[ix] = rng.permutation(x[ix])
            null[i] = stats.spearmanr(resid, xp)[0]
        print(f"   permute within {strat:6s} | {k:8s} rho {obs:+.3f} p {(np.abs(null) >= abs(obs)).mean():.4f}")

# B. recurrence within classes / phyla
def dl(b, se):
    w = 1 / se ** 2; fe = (w * b).sum() / w.sum(); Q = (w * (b - fe) ** 2).sum(); k_ = len(b) - 1
    tau2 = max(0.0, (Q - k_) / (w.sum() - (w ** 2).sum() / w.sum())) if k_ > 0 else 0.0
    wr = 1 / (se ** 2 + tau2); mu = (wr * b).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), max(0.0, (Q - k_) / Q) if Q > 0 else 0.0, stats.chi2.sf(Q, k_) if k_ > 0 else np.nan


for rank, min_n in (("phylum", 40), ("class", 30)):
    out = []
    for grp, gi in B.groupby(rank):
        if len(gi) < min_n:
            continue
        m = (B[rank] == grp).to_numpy()
        try:
            r_, _, _ = fit(m)
        except Exception:
            continue
        out.append({rank: grp, "n": r_["_n"], **{f"{k}_b": r_[k][0] for k in r_ if not k.startswith("_")}, **{f"{k}_se": r_[k][1] for k in r_ if not k.startswith("_")}, **{f"{k}_p": r_[k][2] for k in r_ if not k.startswith("_")}})
    R = pd.DataFrame(out); R.to_csv(D / f"recurrence_by_{rank}_broad.csv", index=False)
    print(f"\nB. recurrence across {len(R)} {rank} groups (n>={min_n})")
    print(R[[rank, "n"] + [c for c in R.columns if c.endswith("_b")]].round(3).to_string(index=False))
    for k, sign in (("Kdp", -1), ("Trk/Ktr", +1), ("His", -1)):
        if f"{k}_b" not in R:
            continue
        ok = R[f"{k}_se"].notna() & R[f"{k}_b"].notna(); b, se = R.loc[ok, f"{k}_b"].to_numpy(), R.loc[ok, f"{k}_se"].to_numpy()
        mu, sem, p, i2, pq = dl(b, se); agree = int((np.sign(b) == sign).sum())
        print(f"   {k}: pooled {mu:+.3f} ± {sem:.3f} (p {p:.1e}), I² {i2:.0%} (Q p {pq:.2f}), expected sign in {agree}/{len(b)} groups (sign test p {stats.binomtest(agree, len(b), .5, alternative='greater').pvalue:.3f})")
