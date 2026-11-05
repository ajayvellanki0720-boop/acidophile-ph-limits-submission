"""Stage 4: model comparison with relatedness-blocked validation.

Outcome: ph_min. Baseline: ph_opt only. Every other model is judged by whether it beats the baseline on
held-out lineages, never by in-sample fit. Random K-fold is deliberately not offered.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

OUTCOME = "ph_min"
COMPOSITION = ["acidic_excess", "mean_pI", "median_pI", "frac_pI_gt7", "frac_pI_lt5", "mean_gravy",
               "frac_Cys", "frac_His", "frac_membrane_proxy", "mean_pI_membrane", "mean_pI_soluble",
               "gc", "genome_size_bp", "ivywrel"]
PROTEOME_CHARGE = ["acidic_excess", "mean_pI", "frac_pI_gt7"]
TEMPERATURE = ["temp_opt_c"]


def design_columns(df: pd.DataFrame, source: str = "hmm"):
    """Resolve feature groups present in df. source='hmm' uses Pfam-validated hm_rate_*/hm_pw_*; 'keyword' uses
    the provisional rate_*/pw_* columns (kept for comparison; detection validation showed them unreliable)."""
    gp, pp = ("hm_rate_", "hm_pw_") if source == "hmm" else ("rate_", "pw_")
    gene = [c for c in df.columns if c.startswith(gp)]
    pathway = [c for c in df.columns if c.startswith(pp)]
    comp = [c for c in COMPOSITION if c in df.columns]
    comp += [c for c in df.columns if c.startswith("rho_")]
    charge = [c for c in PROTEOME_CHARGE if c in df.columns]
    return {"gene": gene, "pathway": pathway, "composition": comp, "charge": charge}


def add_controls(df: pd.DataFrame):
    """One-hot domain and acid identity; return (df, control column names)."""
    out = df.copy()
    ctrl = [c for c in TEMPERATURE if c in out.columns]
    for col in ("domain", "acid_identity"):
        if col in out.columns:
            d = pd.get_dummies(out[col], prefix=col, dtype=float)
            out = pd.concat([out, d], axis=1)
            ctrl += list(d.columns)
    return out, ctrl


def model_specs(df: pd.DataFrame, controls: list, source: str = "hmm"):
    g = design_columns(df, source)
    opt = ["ph_opt"]
    specs = {
        "M1_opt_only": opt,
        "M2_acidic_proteome": g["charge"] + opt,
        "M3_composition": g["composition"] + opt,
        "M4_gene_content": g["gene"] + opt,
        "M5_pathway": g["pathway"] + opt,
        "M6_combined": opt + g["composition"] + g["pathway"] + g["gene"] + controls,
    }
    # Same models without the optimum: asks whether the genome can stand in for ph_opt at all.
    specs.update({k + "_noopt": [c for c in v if c != "ph_opt"]
                  for k, v in list(specs.items()) if k != "M1_opt_only"})
    return {k: v for k, v in specs.items() if v}


def make_splitter(groups: pd.Series, scheme: str, n_splits=5):
    ng = groups.nunique()
    if scheme == "family_blocked":
        return GroupKFold(n_splits=min(n_splits, ng))
    if scheme == "leave_one_lineage_out":
        return LeaveOneGroupOut()
    raise ValueError("scheme must be family_blocked or leave_one_lineage_out")


def _pipeline():
    return make_pipeline(SimpleImputer(strategy="median"), VarianceThreshold(0.0), StandardScaler(),
                         RidgeCV(alphas=np.logspace(-2, 3, 20)))


def cross_val_predict_blocked(df, features, group_col, scheme):
    groups = df[group_col]
    pred = pd.Series(np.nan, index=df.index)
    for tr, te in make_splitter(groups, scheme).split(df, groups=groups):
        assert not set(groups.iloc[tr]) & set(groups.iloc[te]), "relatedness leak across folds"
        m = _pipeline().fit(df.iloc[tr][features], df.iloc[tr][OUTCOME])
        pred.iloc[te] = m.predict(df.iloc[te][features])
    return pred


def _metrics(y, p):
    rmse = float(np.sqrt(np.mean((y - p) ** 2)))
    ss = float(((y - y.mean()) ** 2).sum())
    return {"n": len(y), "rmse": rmse, "mae": float(np.mean(np.abs(y - p))),
            "r2": 1 - float(((y - p) ** 2).sum()) / ss if ss else float("nan"),
            "spearman": float(spearmanr(y, p)[0])}


def group_bootstrap_delta(y, p_base, p_model, groups, n_boot=1000, seed=0):
    """Delta RMSE (model - baseline) with CI from resampling whole groups. Negative = model better."""
    rng = np.random.default_rng(seed)
    gs = pd.Series(groups.values)
    idx = {g: np.where(gs.values == g)[0] for g in gs.unique()}
    keys = list(idx)
    se_b, se_m = (y - p_base).values ** 2, (y - p_model).values ** 2
    deltas = []
    for _ in range(n_boot):
        pick = np.concatenate([idx[k] for k in rng.choice(keys, len(keys))])
        deltas.append(np.sqrt(se_m[pick].mean()) - np.sqrt(se_b[pick].mean()))
    return float(np.sqrt(se_m.mean()) - np.sqrt(se_b.mean())), *np.percentile(deltas, [2.5, 97.5])


def compare_models(df, group_col="family", scheme="family_blocked", min_n=12, source="hmm"):
    """Fit all specs under blocked CV; return metrics table with delta-RMSE vs M1 and its group-bootstrap CI."""
    d, ctrl = add_controls(df.dropna(subset=[OUTCOME, "ph_opt", group_col]).reset_index(drop=True))
    if len(d) < min_n or d[group_col].nunique() < 3:
        raise ValueError(f"need >= {min_n} rows and >= 3 groups; have {len(d)} rows, {d[group_col].nunique()} groups")
    preds = {k: cross_val_predict_blocked(d, f, group_col, scheme) for k, f in model_specs(d, ctrl, source).items()}
    y, base = d[OUTCOME], preds["M1_opt_only"]
    rows = []
    for k, p in preds.items():
        ok = p.notna()
        row = {"model": k, **_metrics(y[ok], p[ok])}
        if k != "M1_opt_only":
            row["delta_rmse"], row["delta_lo"], row["delta_hi"] = group_bootstrap_delta(y[ok], base[ok], p[ok], d.loc[ok, group_col])
            row["beats_baseline"] = row["delta_hi"] < 0
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(preds).assign(strain_id=d["strain_id"], y=y)


def censored_violation_rate(pred: pd.Series, reported_lowest_tested: pd.Series):
    """For 'lowest tested' rows true pH_min <= reported. A prediction above it contradicts the data."""
    return float((pred > reported_lowest_tested).mean())


def compare_by_stratum(df, stratum_col, **kw):
    out = {}
    for s, g in df.groupby(stratum_col):
        try:
            out[s] = compare_models(g, **kw)[0]
        except ValueError as e:
            out[s] = str(e)
    return out
