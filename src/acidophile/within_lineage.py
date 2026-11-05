"""Single-lineage test: does any genomic feature explain pH_min beyond pH_opt among close relatives?

Design (fixed before looking at results):
  - y = pH_min residual after OLS on pH_opt within the lineage; features come from FEATURES (pre-specified list)
  - per-feature statistic: Spearman rho(residual, feature); also leave-one-out change in RMSE when the feature is added
  - null by permuting the residual vector against the genome table (whole genomes keep their feature vectors), 5000 times
  - family-wise control: compare the OBSERVED MAXIMUM |rho| (and best LOO gain) with the permutation distribution of the maximum
  - report detectable effect size: the |rho| a single pre-specified feature would need at alpha=0.05
Lineage relatives are not phylogenetically independent; this test controls the optimum, not tree structure."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata

FEATURES = {
    "charge": ["mean_pI", "acidic_excess", "frac_His"],
    "k_uptake": ["hm_rate_kdp", "hm_rate_trk_ktr"],
    "acid_consumption": ["hm_rate_gad", "hm_rate_adi_cad", "hm_rate_urease"],
    "membrane": ["hm_rate_hopanoid_shc", "hm_rate_cfa"],
    "proton_handling": ["hm_rate_f_atp_synthase", "hm_rate_v_a_atp_synthase", "hm_rate_p_type_atpase"],
    "metal": ["hm_rate_cdf_metal_efflux", "hm_rate_arsenic_resistance"],
    "size": ["genome_size_bp"],
    "thermal_control": ["ivywrel"],
}
ALL_FEATURES = [f for v in FEATURES.values() for f in v]


def _ols_resid(y, x):
    x = np.asarray(x, float)
    X = np.column_stack([np.ones(len(x)), x.reshape(len(x), -1)])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ b, b


def _spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def loo_rmse(y, X):
    """Exact leave-one-out RMSE of OLS with intercept, via the hat matrix."""
    A = np.column_stack([np.ones(len(y)), X])
    H = A @ np.linalg.pinv(A)
    e = y - H @ y
    lev = np.diag(H)
    if (lev >= 1 - 1e-9).any():       # a point fits itself perfectly: LOO undefined, treat as no gain
        return float("inf")
    return float(np.sqrt(np.mean((e / (1 - lev)) ** 2)))


def analyse(df: pd.DataFrame, n_perm=5000, seed=0, features=None, covariates=("ph_opt",)):
    """df needs ph_min, ph_opt and feature columns. Returns (per_feature table, summary dict)."""
    rng = np.random.default_rng(seed)
    feats = [f for f in (features or ALL_FEATURES) if f in df.columns and df[f].nunique() > 1]
    y = df["ph_min"].to_numpy(float)
    cov = df[list(covariates)].to_numpy(float)          # ph_opt, optionally also growth temperature
    opt = cov[:, 0]
    res, beta = _ols_resid(y, cov)
    n = len(df)
    base_loo = loo_rmse(y, cov)
    F = {f: df[f].to_numpy(float) for f in feats}
    zF = {f: (F[f] - F[f].mean()) / F[f].std() for f in feats}

    def stats(r):
        rho = np.array([_spearman(r, F[f]) for f in feats])
        return rho

    obs = stats(res)
    # permutation null: shuffle which genome each residual belongs to
    null_max = np.empty(n_perm)
    null_rho = np.empty((n_perm, len(feats)))
    for i in range(n_perm):
        null_rho[i] = stats(rng.permutation(res))
        null_max[i] = np.abs(null_rho[i]).max()
    p_unc = np.array([(np.abs(null_rho[:, j]) >= abs(obs[j]) - 1e-12).mean() for j in range(len(feats))])
    p_fwer = np.array([(null_max >= abs(o) - 1e-12).mean() for o in obs])
    # LOO gain from adding each feature (OLS opt + feature)
    gain = np.array([loo_rmse(y, np.column_stack([cov, zF[f]])) - base_loo for f in feats])
    tab = pd.DataFrame({"feature": feats, "rho_with_residual": obs, "p_uncorrected": p_unc, "p_family_wise": p_fwer, "delta_loo_rmse": gain})
    tab["group"] = tab.feature.map({f: g for g, v in FEATURES.items() for f in v})
    summary = {"covariates": "+".join(covariates), "n": n, "n_features_tested": len(feats), "opt_slope": float(beta[1]), "opt_r2": float(1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum()),
               "baseline_loo_rmse": base_loo, "residual_sd": float(res.std(ddof=1 + cov.shape[1])),
               "detectable_abs_rho_single": float(np.quantile(np.abs(null_rho[:, 0]), 0.95)),
               "detectable_abs_rho_family_wise": float(np.quantile(null_max, 0.95)),
               "max_abs_rho": float(np.abs(obs).max()), "p_family_wise_max": float((null_max >= np.abs(obs).max() - 1e-12).mean())}
    return tab.sort_values("p_uncorrected").reset_index(drop=True), summary
