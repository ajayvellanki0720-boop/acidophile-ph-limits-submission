"""Broadened-sampling analysis of the histidine and Kdp associations using ALL bacteria with a growth-pH range, with or without a reported optimum. -> DATA_DIR/broadened_by_lineage.csv and DATA_DIR/broadened_contrasts.csv and DATA_DIR/broadened_leave_one_class_out.csv and DATA_DIR/broadened_screen_total_effect.csv
Model (total effect): pH_min ~ source indicators + one standardised feature, PGLS on the marker tree. The optimum is NOT a covariate here (many added species have none);
it is a mediator (about a quarter of the histidine effect ran through it), so this is the total association. A with-optimum version on the subset that has one is run for comparison.
 A. screen of all features (Benjamini-Hochberg) on everything
 B. per-phylum and per-class histidine / Kdp estimates on the comparable scale (pH units per percentage point of histidine; per 0.1 Kdp copies per 1000 genes) + pooling
 C. independent contrasts per phylum (all contrasts)
 D. leave-one-class-out confirmation
Usage: python3 scripts/29_total_effect_broadened.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
from acidophile.within_lineage import ALL_FEATURES

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
B["phylum"] = B.phylum.fillna("NA"); B["class"] = B["class"].fillna("NA")
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
FEATS = [f for f in ALL_FEATURES + ["gc"] if f in B and B[f].std() > 0]
C_ALL = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float); Fm = B[FEATS].to_numpy(float)
print(f"bacteria: {len(B)} species | with optimum {B.ph_opt.notna().sum()} | pH_min {B.ph_min.min():.1f}-{B.ph_min.max():.1f}, SD {B.ph_min.std():.2f} | phyla {B.phylum.nunique()} | origin {B.origin.value_counts().to_dict()}")


def covs(g, with_opt=False):
    return (["ph_opt"] if with_opt else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]


def scan(idx, feats=None, with_opt=False):
    g = B.iloc[idx]
    X = np.column_stack([np.ones(len(idx))] + [g[c].to_numpy(float) for c in covs(g, with_opt)])
    cols = [FEATS.index(f) for f in (feats or FEATS)]
    b, se, p, lam = ph.pgls_scan(C_ALL[np.ix_(idx, idx)], Y[idx], X, Fm[:, cols][idx])
    return pd.DataFrame({"feature": [FEATS[c] for c in cols], "beta": b, "se": se, "p": p, "n": len(idx), "lambda": lam})


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]; return np.minimum(q, 1)


def dl(b, se):
    w = 1 / se ** 2; fe = (w * b).sum() / w.sum(); Q = (w * (b - fe) ** 2).sum(); k = len(b) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (se ** 2 + tau2); mu = (wr * b).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), (max(0.0, (Q - k) / Q) if Q > 0 else 0.0), (stats.chi2.sf(Q, k) if k > 0 else np.nan)


allidx = np.arange(len(B))
A = scan(allidx); A["q_BH"] = bh(A.p.fillna(1)); A = A.sort_values("p"); A.to_csv(D / "broadened_screen_total_effect.csv", index=False)
print("\nA. screen, total effect, all bacteria (beta = pH_min units per SD of the feature):"); print(A.head(8).round(5).to_string(index=False))
wo = np.where(B.ph_opt.notna())[0]
A2 = scan(wo, with_opt=True); A2["q_BH"] = bh(A2.p.fillna(1)); print("\n   same screen on the subset with an optimum, optimum as covariate (n=%d):" % len(wo)); print(A2.sort_values("p").head(5).round(5).to_string(index=False))

print("\nB. per-lineage estimates on the comparable scale")
rows = []
for rank, min_n in (("phylum", 40), ("class", 40)):
    for grp, gi in B.groupby(rank):
        if len(gi) < min_n or grp == "NA":
            continue
        idx = np.where(B[rank] == grp)[0]
        r = scan(idx, ["frac_His", "hm_rate_kdp"]).set_index("feature")
        sd_h, sd_k = B.iloc[idx].frac_His.std() * 100, B.iloc[idx].hm_rate_kdp.std()
        rows.append({"rank": rank, "lineage": grp, "n": len(idx), "His_SD_pts": sd_h, "His_pH_per_pt": r.loc["frac_His", "beta"] / sd_h, "His_se": r.loc["frac_His", "se"] / sd_h, "His_p": r.loc["frac_His", "p"],
                     "Kdp_pH_per_0.1": 0.1 * r.loc["hm_rate_kdp", "beta"] / sd_k, "Kdp_se": 0.1 * r.loc["hm_rate_kdp", "se"] / sd_k, "Kdp_p": r.loc["hm_rate_kdp", "p"]})
L = pd.DataFrame(rows); L.to_csv(D / "broadened_by_lineage.csv", index=False)
for rank in ("phylum", "class"):
    x = L[L["rank"] == rank]; print(f"\n  {rank} (n>=40):"); print(x.drop(columns="rank").round(3).to_string(index=False))
    for lab, bcol, scol in (("His, pH per +1 pt", "His_pH_per_pt", "His_se"), ("Kdp, pH per +0.1 rate", "Kdp_pH_per_0.1", "Kdp_se")):
        mu, sem, p, i2, pq = dl(x[bcol].to_numpy(), x[scol].to_numpy()); print(f"   pooled {lab}: {mu:+.3f} ± {sem:.3f} (p {p:.1e}); I² {i2:.0%} (Q p {pq:.2f}); negative in {(x[bcol] < 0).sum()}/{len(x)} lineages")

print("\nC. independent contrasts per phylum (all contrasts; per pH unit per +1 pt of histidine)")
rows = []
for grp in ("Bacillota", "Pseudomonadota", "Actinomycetota", "Bacteroidota"):
    g = B[B.phylum == grp].reset_index(drop=True)
    cols = [c for c in ("src_thermobase", "src_acidb", "src_bacdive") if g[c].nunique() > 1]
    sd_h = g.frac_His.std() * 100
    M = np.column_stack([g.ph_min] + [g[c] for c in cols] + [(g.frac_His - g.frac_His.mean()) / g.frac_His.std(), (g.hm_rate_kdp - g.hm_rate_kdp.mean()) / g.hm_rate_kdp.std()])
    c, ch = ph.pic_contrasts(tree, list(g.ncbi_assembly), M)
    b, se, p = ph.regress_through_origin(c[:, 0], c[:, 1:]); k = len(cols)
    rows.append({"phylum": grp, "n": len(g), "His_SD_pts": sd_h, "His_pH_per_pt": b[k] / sd_h, "se": se[k] / sd_h, "p": p[k], "Kdp_beta_per_sd": b[k + 1], "Kdp_se": se[k + 1], "Kdp_p": p[k + 1]})
Cc = pd.DataFrame(rows); Cc.to_csv(D / "broadened_contrasts.csv", index=False); print(Cc.round(3).to_string(index=False))
mu, sem, p, i2, pq = dl(Cc.His_pH_per_pt.to_numpy(), Cc.se.to_numpy()); print(f"   pooled His {mu:+.3f} ± {sem:.3f} pH per pt (p {p:.3f}); I² {i2:.0%} (Q p {pq:.2f})")

print("\nD. leave-one-class-out confirmation (training screen on all other bacteria, BH q<0.05; test in the held-out class)")
rows = []
for grp in [k for k, v in B.groupby("class").size().items() if v >= 40 and k != "NA"]:
    te = np.where(B["class"] == grp)[0]; tr = np.where(B["class"] != grp)[0]
    S = scan(tr); S["q"] = bh(S.p.fillna(1)); T_ = scan(te)
    m = S.merge(T_, on="feature", suffixes=("_tr", "_te")); m["fold"] = grp; rows.append(m)
M = pd.concat(rows); M.to_csv(D / "broadened_leave_one_class_out.csv", index=False)
for f in M[M.q < 0.05].feature.unique():
    x = M[M.feature == f]; sel = x[x.q < 0.05]; sign = np.sign(sel.beta_tr.median())
    mu, sem, p, i2, _ = dl(sel.beta_te.to_numpy(), sel.se_te.to_numpy()); agree = int((np.sign(sel.beta_te) == sign).sum())
    print(f"  {f:22s} selected in {len(sel)}/{len(x)} folds | held-out pooled {mu:+.3f} ± {sem:.3f} (p {p:.1e}, I² {i2:.0%}) | same sign {agree}/{len(sel)}")
