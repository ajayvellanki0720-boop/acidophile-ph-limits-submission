"""Lineage hold-out confirmation test. Discovery (feature screen with Benjamini-Hochberg) and confirmation use DIFFERENT lineages, so the
confirmation estimate is not affected by selection. Bacteria only (Kdp, Trk/Ktr and His are bacterial leads; archaea tested separately at the end).
 1. leave-one-class-out   2. leave-one-phylum-out   3. random half-splits of classes (replication rate per feature, all features run identically)
Model: PGLS pH_min ~ pH_opt + source indicators (those varying in the subset) + one standardised feature. In each subset Pagel's lambda is fitted at the
covariates-only baseline and then held fixed for the feature screen (closed-form, via the whitened design), which is what makes hundreds of splits feasible.
Usage: python3 scripts/26_lineage_holdout.py [n_random_splits]   (env DATA_DIR, PHYLO_DIR)"""
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
from scipy.linalg import solve_triangular
from acidophile import phylo as ph
from acidophile.within_lineage import ALL_FEATURES

warnings.filterwarnings("ignore")
NSPLIT = int(sys.argv[1]) if len(sys.argv) > 1 else 300
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
B["class"] = B["class"].fillna("NA"); B["phylum"] = B["phylum"].fillna("NA")
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    B[f"src_{o}"] = (B.origin == o).astype(float)
FEATS = [f for f in ALL_FEATURES + ["gc"] if f in B and B[f].std() > 0]
CALL = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float)
F = B[FEATS].to_numpy(float)
rng = np.random.default_rng(2026)


def screen(idx, want=None):
    """Return DataFrame(feature, beta, se, p) for the subset `idx` (array of row indices); lambda fixed at the covariates-only baseline."""
    g = B.iloc[idx]
    cov = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive", "src_literature_auto") if g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(idx))] + [g[c].to_numpy(float) for c in cov]); C = CALL[np.ix_(idx, idx)]
    base = ph.pgls(Y[idx], X, C); lam = base["lambda"]
    s = np.mean(np.diag(C)); Cn = C / s
    V = Cn * lam + np.diag(np.diag(Cn)) * (1 - lam)
    L = np.linalg.cholesky(V + 1e-10 * np.eye(len(idx)))
    Xs = solve_triangular(L, X, lower=True); ys = solve_triangular(L, Y[idx], lower=True)
    Q, _ = np.linalg.qr(Xs)
    yt = ys - Q @ (Q.T @ ys)
    Fz = F[idx]; sd = Fz.std(axis=0); ok = sd > 0
    Fz = (Fz - Fz.mean(axis=0)) / np.where(ok, sd, 1)
    Fs = solve_triangular(L, Fz, lower=True); Ft = Fs - Q @ (Q.T @ Fs)
    xx = (Ft ** 2).sum(axis=0); xx[~ok] = np.nan
    beta = (Ft * yt[:, None]).sum(axis=0) / xx
    dfree = len(idx) - Xs.shape[1] - 1
    s2 = ((yt ** 2).sum() - beta ** 2 * xx) / dfree
    se = np.sqrt(s2 / xx); p = 2 * stats.t.sf(np.abs(beta / se), dfree)
    return pd.DataFrame({"feature": FEATS, "beta": beta, "se": se, "p": p, "n": len(idx), "lambda": lam})


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty_like(p); q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]; return np.minimum(q, 1)


def dl(b, se):
    w = 1 / se ** 2; fe = (w * b).sum() / w.sum(); Q = (w * (b - fe) ** 2).sum(); k = len(b) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (se ** 2 + tau2); mu = (wr * b).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), (max(0.0, (Q - k) / Q) if Q > 0 else 0.0)


FOCUS = ["frac_His", "hm_rate_kdp", "hm_rate_trk_ktr"]
out_rows, fold_rows = [], []
for rank, min_n in (("class", 25), ("phylum", 25)):
    groups = [k for k, v in B.groupby(rank).size().items() if v >= min_n and k != "NA"]
    folds = []
    for gname in groups:
        te = np.where(B[rank] == gname)[0]; tr = np.where(B[rank] != gname)[0]
        S_tr = screen(tr); S_tr["q"] = bh(S_tr.p.fillna(1)); S_te = screen(te)
        m = S_tr.merge(S_te, on="feature", suffixes=("_train", "_test"))
        m["fold"] = gname; m["rank"] = rank; folds.append(m)
    M = pd.concat(folds); fold_rows.append(M)
    print(f"\n=== leave-one-{rank}-out: {len(groups)} held-out groups (n>={min_n}) ===")
    for f in FEATS:
        x = M[M.feature == f]; sel = x[x.q < 0.05]
        if len(sel) == 0:
            continue
        sign = np.sign(sel.beta_train.median()); b, se = sel.beta_test.to_numpy(), sel.se_test.to_numpy()
        mu, sem, pp, i2 = dl(b, se); agree = int((np.sign(b) == sign).sum())
        z = (b * sign / se); stouffer = stats.norm.sf(z.sum() / np.sqrt(len(z)))
        out_rows.append({"scheme": f"leave-one-{rank}-out", "feature": f, "folds_selected": f"{len(sel)}/{len(x)}", "train_sign": int(sign), "heldout_pooled_beta": mu, "se": sem, "p_two_sided": pp,
                         "I2": i2, "heldout_same_sign": f"{agree}/{len(sel)}", "one_sided_stouffer_p": stouffer})
        print(f"  {f:26s} selected in {len(sel):2d}/{len(x)} training folds (sign {int(sign):+d}) | held-out pooled {mu:+.3f} ± {sem:.3f} (p {pp:.1e}, I² {i2:.0%}) | same sign {agree}/{len(sel)} | one-sided Stouffer p {stouffer:.2e}")
pd.concat(fold_rows).to_csv(D / "holdout_leave_one_lineage_out_folds.csv", index=False)

# 3. random half-splits of classes
cls = B.groupby("class").size(); cls = cls[cls.index != "NA"]; names = np.array(cls.index)
rec = []
done = 0
while done < NSPLIT:
    perm = rng.permutation(len(names)); half = np.zeros(len(names), bool); half[perm[: len(names) // 2]] = True
    A = np.where(B["class"].isin(names[half]))[0]; Bi = np.where(B["class"].isin(names[~half]) | (B["class"] == "NA"))[0]
    if min(len(A), len(Bi)) < 450 or max(len(A), len(Bi)) > 850:
        continue                                        # balanced halves only
    SA = screen(A); SA["q"] = bh(SA.p.fillna(1)); SB = screen(Bi)
    m = SA.merge(SB, on="feature", suffixes=("_disc", "_conf")); m["split"] = done; rec.append(m); done += 1
R = pd.concat(rec); R.to_csv(D / "holdout_random_class_splits.csv", index=False)
rows = []
for f in FEATS:
    x = R[R.feature == f]; d = x[x.q < 0.05]
    if len(d) == 0:
        rows.append({"feature": f, "discovered": 0, "disc_rate": 0.0}); continue
    sign = np.sign(d.beta_disc); zc = d.beta_conf * sign / d.se_conf
    rows.append({"feature": f, "discovered": len(d), "disc_rate": len(d) / len(x), "replicated_p05_onesided": int((zc > 1.645).sum()), "replication_rate": float((zc > 1.645).mean()),
                 "same_sign_rate": float((np.sign(d.beta_conf) == sign).mean()), "median_conf_beta_signed": float((d.beta_conf * sign).median())})
S = pd.DataFrame(rows).sort_values("disc_rate", ascending=False); S.to_csv(D / "holdout_random_split_summary.csv", index=False)
print(f"\n=== random half-splits of classes ({NSPLIT} splits, discovery q<0.05 in one half, one-sided test in the other) ===")
print(S[S.discovered > 0].round(3).to_string(index=False))
pd.DataFrame(out_rows).to_csv(D / "holdout_leave_one_lineage_out_summary.csv", index=False)
