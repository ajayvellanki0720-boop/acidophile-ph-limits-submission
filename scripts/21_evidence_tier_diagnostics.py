"""Does the genome-based signal depend on the evidence tier of the minima? Reproduces, in one place:
 (1) PGLS Kdp / Trk-Ktr by tier on the bacterial ML tree, (2) K+ pair added to the optimum under family-blocked CV by tier,
 (3) genome-only R2 on random subsets of the unverified rows at the verified sample size, (4) temperature + domain as competing predictors,
 (5) tier structure (optimum spread) and gene content -> optimum. Usage: python3 scripts/21_evidence_tier_diagnostics.py
Tiers: verified = literature species description with a type-matched genome (confidence medium/high, all literature, no data-source mixing);
unverified = ThermoBase / AciDB / BacDive rows and low-confidence literature rows."""
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
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/full")      # DATA_DIR=results/broad for the expanded set
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[df.ph_opt.notna()].reset_index(drop=True)
df["verified"] = df.confidence.isin(["medium", "high"]); df["archaea"] = (df.domain == "Archaea").astype(float); df["temp_known"] = df.temp_opt_c.notna().astype(float)
gene = [c for c in df.columns if c.startswith("hm_rate_") and df[c].nunique() > 1]
z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)


def blocked_pred(sub, cols, target, seed):
    fam = sub.family.fillna("NA"); order = np.random.default_rng(seed).permutation(fam.unique()); g = fam.map({f: i % 8 for i, f in enumerate(order)}).to_numpy(); p = np.full(len(sub), np.nan)
    for k in range(8):
        tr, te = np.where(g != k)[0], np.where(g == k)[0]
        if len(te) == 0 or len(tr) < 10:
            continue
        m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), RidgeCV(alphas=np.logspace(-1, 4, 20))).fit(sub.loc[tr, cols], sub.loc[tr, target]); p[te] = m.predict(sub.loc[te, cols])
    return p


def r2(sub, cols, target="ph_min", seeds=4):
    sub = sub.reset_index(drop=True); y = sub[target].to_numpy(); out = []
    for s in range(seeds):
        p = blocked_pred(sub, cols, target, s); out.append(1 - np.nansum((y - p) ** 2) / np.sum((y - y.mean()) ** 2))
    return float(np.mean(out))


V, U = df[df.verified], df[~df.verified]
# (1) PGLS by tier (bacteria, ML tree)
tree = Phylo.read(ROOT / os.environ.get("PHYLO_DIR", "results/phylo_ml") / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
rows = []
for label, mask in (("verified", B.verified), ("unverified", ~B.verified), ("ThermoBase only", B.origin == "thermobase"), ("all bacteria", np.ones(len(B), bool))):
    g = B[mask].reset_index(drop=True)
    if len(g) < 25 or g.hm_rate_kdp.std() == 0 or g.hm_rate_trk_ktr.std() == 0:
        continue
    cov = ["ph_opt"]
    for o in ("thermobase", "acidb", "bacdive"):
        g[f"s_{o}"] = (g.origin == o).astype(float)
        if g[f"s_{o}"].nunique() > 1:
            cov.append(f"s_{o}")
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in cov] + [z(g.hm_rate_kdp), z(g.hm_rate_trk_ktr)])
    r = ph.pgls(g.ph_min.to_numpy(float), X, ph.vcv_fast(tree, list(g.ncbi_assembly))); k = X.shape[1] - 2
    rows.append({"tier": label, "n": len(g), "kdp_b": r["beta"][k], "kdp_se": r["se"][k], "kdp_p": r["p"][k], "trk_b": r["beta"][k + 1], "trk_se": r["se"][k + 1], "trk_p": r["p"][k + 1]})
P = pd.DataFrame(rows); P.to_csv(D / "kplus_by_evidence_tier.csv", index=False)
v, u = P.iloc[0], P.iloc[1]
print("(1) joint PGLS by tier:\n", P.round(4).to_string(index=False))
for f in ("kdp", "trk"):
    dif = u[f + "_b"] - v[f + "_b"]; se = np.hypot(u[f + "_se"], v[f + "_se"]); print(f"    {f}: unverified minus verified {dif:+.3f} (z {dif / se:+.2f}, p {2 * stats.norm.sf(abs(dif / se)):.2f})")
# (3) subsample test
rv = r2(V, gene)
draws = [r2(U.sample(len(V), random_state=i), gene, seeds=2) for i in range(40)]
print(f"(3) gene-content-only R2: verified n={len(V)} {rv:.2f}; {len(V)}-species draws from unverified: median {np.median(draws):.2f}, 5-95% {np.percentile(draws, 5):.2f} to {np.percentile(draws, 95):.2f}")
# (4) temperature + domain
T = []
for name, sub in (("all", df), ("verified", V), ("unverified", U)):
    T.append({"subset": name, "n": len(sub), "temperature+domain": r2(sub, ["temp_opt_c", "archaea", "temp_known"]), "gene content": r2(sub, gene), "optimum only": r2(sub, ["ph_opt"])})
T = pd.DataFrame(T).set_index("subset"); T.to_csv(D / "genome_vs_temperature_domain_r2.csv"); print("(4) family-blocked R2 for pH_min:\n", T.round(2).to_string())
# (5) structure
S = pd.DataFrame([{"tier": n, "n": len(s), "pH_opt sd": s.ph_opt.std(), "pH_min sd": s.ph_min.std(), "mean temp_opt": s.temp_opt_c.mean(), "archaea share": (s.domain == "Archaea").mean(),
                   "gene->pH_opt R2": r2(s, gene, "ph_opt"), "gene->pH_min R2": r2(s, gene, "ph_min")} for n, s in (("verified", V), ("unverified", U))]).set_index("tier")
S.to_csv(D / "tier_structure_and_genome_prediction.csv"); print("(5)\n", S.round(2).T.to_string())
# (2) K+ pair in blocked CV, bacteria
Bd = df[df.domain == "Bacteria"].reset_index(drop=True)
rng = np.random.default_rng(5)
out2 = []
for label, sub in (("verified", Bd[Bd.verified]), ("unverified", Bd[~Bd.verified]), ("all", Bd)):
    sub = sub.reset_index(drop=True); fam = sub.family.fillna("NA"); idx = {f: np.where(fam == f)[0] for f in fam.unique()}; keys = list(idx); res = []
    for s in range(5):
        e0 = (sub.ph_min - blocked_pred(sub, ["ph_opt"], "ph_min", s)) ** 2; e1 = (sub.ph_min - blocked_pred(sub, ["ph_opt", "hm_rate_kdp", "hm_rate_trk_ktr"], "ph_min", s)) ** 2
        dd = [np.sqrt(e1.values[np.concatenate([idx[k] for k in rng.choice(keys, len(keys))])].mean()) - np.sqrt(e0.values[np.concatenate([idx[k] for k in rng.choice(keys, len(keys))])].mean()) for _ in range(0)]
        pick = [np.concatenate([idx[k] for k in rng.choice(keys, len(keys))]) for _ in range(1000)]
        dd = [np.sqrt(e1.values[p].mean()) - np.sqrt(e0.values[p].mean()) for p in pick]
        res.append((np.sqrt(e1.mean()) - np.sqrt(e0.mean()), np.percentile(dd, 2.5), np.percentile(dd, 97.5)))
    r = np.array(res); out2.append({"tier": label, "n": len(sub), "delta_rmse_mean": r[:, 0].mean(), "ci_lo_median": np.median(r[:, 1]), "ci_hi_median": np.median(r[:, 2]), "seeds_ci_excludes_0": int((r[:, 2] < 0).sum())})
Q = pd.DataFrame(out2); Q.to_csv(D / "kplus_cv_gain_by_tier.csv", index=False); print("(2) K+ pair added to optimum (bacteria):\n", Q.round(3).to_string(index=False))
