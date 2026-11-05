"""Is the change point stable enough for the design document to rest on? -> DATA_DIR/changepoint_robustness.csv
scripts/63 found a hinge in the composite that beats both the line and the quadratic on AIC, put the break at +0.217 SD, and concluded that the
proposed experiment's chassis sits below it and therefore inside the live part of the range. docs/CHARGE_REVERSAL_DESIGN.md section 4c now rests on
that conclusion, and the break was fitted once, on one composite, with one set of covariates. A break point is a notoriously unstable thing to
estimate - the likelihood surface is flat over wide stretches when the slope change is modest, and the profile interval from scripts/63 already runs
from -0.147 to +0.548, which is a third of the predictor's range. So the claim needs checking before it is allowed to stand, and the check has to be
aimed at the claim that is actually being made. That claim is not "the break is at +0.217". It is "E. coli is on the low side of the break", which is a
much weaker and more robust statement, and the right test is how often that survives rather than how tightly the break is pinned.
Two parts.
  1. SPECIFICATION VARIANTS. The break is refitted on the factor-weighted composite, on the 19 high-coverage references, with the growth optimum added
     as a covariate, on the topology-only periplasmic composite from scripts/59, and with the 298 rows that scripts/61 flags as optimum-read-as-range
     dropped. Each refits lambda by maximum likelihood. For each, the break, the two slopes, the AIC against the line, and where the chassis falls.
     The chassis has to be relocated in every variant, because each composite standardises against its own references and its own sample.
  2. HALF-SAMPLE STABILITY. A family bootstrap is not available here - resampling tips with replacement makes the phylogenetic covariance singular -
     so the sampling behaviour of the break is estimated by repeatedly splitting the 428 families at random into two halves and refitting the break in
     each half. That is a subsampling estimator rather than a bootstrap, and it understates the spread a full-sample estimate would have, which is
     stated rather than hidden. Lambda is held at the full-sample value inside the splits so that one Cholesky factorisation serves all 41 grid points
     of a split; that is the same cost concession scripts/50 made, and it removes a source of variability, so the spread reported here is again
     optimistic. What matters is the fraction of splits in which the chassis still falls below the fitted break.
Usage: python3 scripts/64_changepoint_robustness.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.linalg import solve_triangular
from acidophile import phylo as ph
from acidophile.parser_flags import flagged_pmids

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
rng = np.random.default_rng(0)
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
cls_tab = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv")
reflen = cls_tab.groupby(["ref", "class"]).size()
wide = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
W = wide.set_index("ncbi_assembly")
_used = set(B[B.origin == "literature_auto"].source.dropna().str.replace("PMID:", "", regex=False).astype(str))
_bad, _ = flagged_pmids(pd.read_csv(ROOT / "data/curated/broad_candidates.csv"), _used)
flagged = set(B.ncbi_assembly[B.source.dropna().str.replace("PMID:", "", regex=False).astype(str).isin(_bad).reindex(B.index, fill_value=False)])
print(f"{len(B)} bacteria on the tree; {len(flagged)} of them carry the optimum-as-range flag (same detector as scripts/61)")


def per_ref(label, classes):
    """Acidic excess over the union of classes for one reference, recombined on counts as scripts/59 does."""
    num = pd.Series(0.0, index=B.ncbi_assembly)
    den = pd.Series(0.0, index=B.ncbi_assembly)
    ok = pd.Series(False, index=B.ncbi_assembly)
    floor = 0.0
    for cls in classes:
        c = f"{label}|{cls}|n"
        if c not in W or (label, cls) not in reflen:
            continue
        num = num + (W[f"{label}|{cls}|D"] + W[f"{label}|{cls}|E"] - W[f"{label}|{cls}|K"] - W[f"{label}|{cls}|R"]).reindex(B.ncbi_assembly).fillna(0.0)
        den = den + W[c].reindex(B.ncbi_assembly).fillna(0.0)
        ok = ok | W[c].reindex(B.ncbi_assembly).notna()
        floor += 0.5 * reflen[(label, cls)]
    if not ok.any():
        return None, np.nan
    x = (num / den).where(ok & (den >= max(8, floor)))
    # the reference's own sequence, which for every wide-panel protein is E. coli's
    s = cls_tab[(cls_tab.ref == label) & cls_tab["class"].isin(classes)]
    own = np.nan
    if len(s) >= 8:
        aa = s.aa.to_numpy()
        own = (np.isin(aa, ["D", "E"]).sum() - np.isin(aa, ["K", "R"]).sum()) / len(aa)
    return pd.Series(x.to_numpy(), index=B.index), own


def composite(classes=("peri_exp",), cover=0.0, weights=None):
    """Returns the composite, the chassis's value on the same scale, and the reference count."""
    cols, chassis = {}, {}
    for r in MEM:
        x, own = per_ref(r, list(classes))
        if x is None or x.notna().mean() < cover:
            continue
        mu, sd = float(x.mean()), float(x.std())
        cols[r] = (x - mu) / sd
        if np.isfinite(own):
            chassis[r] = (own - mu) / sd
    Z = pd.DataFrame(cols)
    if weights == "pc":
        cc = Z.notna().all(axis=1)
        w = np.linalg.eigh(np.corrcoef(Z[cc].to_numpy(float), rowvar=False))[1][:, -1]
        if w.sum() < 0:
            w = -w
        Wt = pd.DataFrame(np.tile(w, (len(Z), 1)), columns=Z.columns, index=Z.index).where(Z.notna())
        s = (Z * Wt).sum(axis=1, min_count=3) / Wt.abs().sum(axis=1)
        wmap = dict(zip(Z.columns, np.abs(w)))
        ch = sum(chassis[r] * wmap[r] for r in chassis) / sum(wmap[r] for r in chassis)
    else:
        s = Z.mean(axis=1)
        ch = float(np.mean([chassis[r] for r in chassis]))
    return s.where(Z.notna().sum(axis=1) >= 3), ch, len(cols)


GRID = np.linspace(0.10, 0.90, 41)


def hinge(y, X0, x, C, lam=None):
    """Profile the break over quantiles of x. Returns (break, b_low, b_extra, se_extra, loglik, best index)."""
    grid = np.quantile(x, GRID)
    out = []
    for c in grid:
        f = ph.pgls(y, np.column_stack([X0, x, np.maximum(0.0, x - c)]), C, lam=lam)
        out.append((c, float(f["loglik"]), float(f["beta"][-2]), float(f["beta"][-1]), float(f["se"][-1])))
    P = pd.DataFrame(out, columns=["c", "loglik", "b_low", "b_extra", "se_extra"])
    return P.loc[P.loglik.idxmax()], P


VARIANTS = [
    ("the published fit (32 refs, exposed periplasmic)", dict(), ["gw"], False),
    ("factor-weighted, 19 high-coverage", dict(cover=0.80, weights="pc"), ["gw"], False),
    ("unweighted, 19 high-coverage", dict(cover=0.80), ["gw"], False),
    ("topology only, periplasmic all positions", dict(classes=("peri_exp", "peri_bur")), ["gw"], False),
    ("growth optimum also controlled", dict(), ["gw", "ph_opt"], False),
    ("optimum-as-range rows dropped", dict(), ["gw"], True),
]
rows = []
print("\n1. specification variants (break in SD of that variant's own composite)")
print(f"   {'variant':46s} {'refs':>4} {'n':>5} {'break':>7} {'below':>8} {'above':>8} {'dAIC vs line':>13} {'chassis':>8} side")
store = {}
for name, kw, covs, drop_flagged in VARIANTS:
    pr, chassis, k = composite(**kw)
    m = pr.notna() & B.ph_min.notna() & B[covs].notna().all(axis=1)
    if drop_flagged:
        m = m & ~B.ncbi_assembly.isin(flagged)
    m = m.to_numpy()
    g = B[m]
    x_raw = pr[m].to_numpy(float)
    mu, sd = x_raw.mean(), x_raw.std()
    x = (x_raw - mu) / sd
    ch = (chassis - mu) / sd
    y = g.ph_min.to_numpy(float)
    use = [c for c in covs + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X0 = np.column_stack([np.ones(len(y))] + [g[c].to_numpy(float) for c in use])
    C = Vcv[np.ix_(m, m)]
    lin = ph.pgls(y, np.column_stack([X0, x]), C)
    best, P = hinge(y, X0, x, C)
    aic_h = -2 * best.loglik + 2 * (X0.shape[1] + 3)
    aic_l = -2 * float(lin["loglik"]) + 2 * (X0.shape[1] + 1)
    within = P[P.loglik >= best.loglik - 1.92]
    side = "BELOW" if ch < best.c else "above"
    rows.append({"part": "variant", "variant": name, "references": k, "n": int(m.sum()), "break_sd": float(best.c),
                 "break_lo": float(within.c.min()), "break_hi": float(within.c.max()), "slope_low": float(best.b_low),
                 "slope_high": float(best.b_low + best.b_extra), "se_extra": float(best.se_extra), "aic_hinge": aic_h,
                 "aic_linear": aic_l, "chassis_per_sd": float(ch), "chassis_side": side.lower(),
                 "chassis_below_interval": bool(ch < within.c.min())})
    print(f"   {name:46s} {k:>4d} {int(m.sum()):>5d} {best.c:>+7.3f} {best.b_low:>+8.4f} {best.b_low + best.b_extra:>+8.4f} "
          f"{aic_h - aic_l:>+13.2f} {ch:>+8.3f} {side}" + ("" if ch >= within.c.min() else ", below its whole interval"))
    if name.startswith("the published"):
        store = dict(y=y, X0=X0, x=x, C=C, fam=g.family.to_numpy(), ch=ch, lam=float(lin["lambda"]))
print("   A negative dAIC means the hinge is preferred. The chassis column is E. coli on that variant's own standardised scale,")
print("   which moves between variants because each composite standardises against its own references and sample.")

print("\n2. half-sample stability of the break, and of the claim that rests on it")
y, X0, x, C, fam, ch, lam = (store[k] for k in ("y", "X0", "x", "C", "fam", "ch", "lam"))
fams = np.unique(fam)
grid = np.quantile(x, GRID)
NSPLIT = 200
breaks, below, sl, sh = [], 0, [], []
for i in range(NSPLIT):
    pick = set(rng.choice(fams, size=len(fams) // 2, replace=False))
    s = np.array([f in pick for f in fam])
    if s.sum() < 400:
        continue
    Cs = C[np.ix_(s, s)]
    V = Cs * lam + np.diag(np.diag(Cs)) * (1 - lam)
    L = np.linalg.cholesky(V + 1e-10 * np.eye(int(s.sum())))
    ys = solve_triangular(L, y[s], lower=True)
    best_ll, best_c, best_b = -np.inf, np.nan, (np.nan, np.nan)
    for c in grid:
        X = np.column_stack([X0[s], x[s], np.maximum(0.0, x[s] - c)])
        Xs = solve_triangular(L, X, lower=True)
        beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
        r = ys - Xs @ beta
        ll = -0.5 * int(s.sum()) * np.log(r @ r)   # concentrated Gaussian log-likelihood, up to constants shared across the grid
        if ll > best_ll:
            best_ll, best_c, best_b = ll, c, (beta[-2], beta[-1])
    breaks.append(best_c)
    sl.append(best_b[0])
    sh.append(best_b[0] + best_b[1])
    below += int(ch < best_c)
breaks = np.asarray(breaks)
print(f"   {len(breaks)} random half-samples of the 428 families, lambda held at the full-sample value {lam:.3f}")
print(f"   break: median {np.median(breaks):+.3f}, 10th to 90th percentile {np.percentile(breaks, 10):+.3f} to {np.percentile(breaks, 90):+.3f}, "
      f"full range {breaks.min():+.3f} to {breaks.max():+.3f}")
print(f"   slope below the break: median {np.median(sl):+.4f};  above: median {np.median(sh):+.4f}")
print(f"   the chassis ({ch:+.3f}) falls below the fitted break in {below}/{len(breaks)} = {below / len(breaks):.1%} of half-samples")
rows.append({"part": "halfsample", "variant": "published fit", "n_splits": len(breaks), "break_median": float(np.median(breaks)),
             "break_p10": float(np.percentile(breaks, 10)), "break_p90": float(np.percentile(breaks, 90)),
             "break_min": float(breaks.min()), "break_max": float(breaks.max()), "slope_low_median": float(np.median(sl)),
             "slope_high_median": float(np.median(sh)), "chassis_per_sd": float(ch), "chassis_below_frac": below / len(breaks)})

R = pd.DataFrame(rows)
R.to_csv(D / "changepoint_robustness.csv", index=False)
print("\nwrote", D / "changepoint_robustness.csv")
print("\nReading. The break's own location is the fragile quantity and should not be quoted to three decimals anywhere. The claim the")
print("design document rests on is the ordinal one - the chassis is on the low side - and the right way to report it is the fraction of")
print("variants and of half-samples in which that holds, not the point estimate of the break.")
