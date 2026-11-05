"""Stage 5: matched close-relative contrasts.

Pairs of relatives with similar pH_opt but different pH_min isolate the machinery that sets the
boundary from general acid adaptation. Tree-based reconstruction of independent increases in acid
tolerance waits for the marker tree; until then pairs are formed within a taxonomic rank.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon


def matched_pairs(df, rank="genus", opt_tol=0.5, min_min_diff=1.0):
    rows = []
    for taxon, g in df.dropna(subset=["ph_opt", "ph_min"]).groupby(rank):
        for a, b in combinations(g.index, 2):
            if abs(g.at[a, "ph_opt"] - g.at[b, "ph_opt"]) <= opt_tol and abs(g.at[a, "ph_min"] - g.at[b, "ph_min"]) >= min_min_diff:
                tol, sens = (a, b) if g.at[a, "ph_min"] < g.at[b, "ph_min"] else (b, a)
                rows.append({"taxon": taxon, "more_tolerant": g.at[tol, "strain_id"], "less_tolerant": g.at[sens, "strain_id"],
                             "i_tol": tol, "i_sens": sens,
                             "d_ph_min": g.at[tol, "ph_min"] - g.at[sens, "ph_min"]})
    return pd.DataFrame(rows)


def contrast_features(df, pairs, features, min_pairs=6):
    """Per feature: mean paired difference (tolerant - sensitive), sign test, Wilcoxon.
    The unit is the taxon: pairs inside one taxon are averaged first so one genus cannot dominate."""
    if pairs.empty:
        return pd.DataFrame()
    out = []
    for f in features:
        d = pairs.assign(diff=[df.at[r.i_tol, f] - df.at[r.i_sens, f] for r in pairs.itertuples()])
        per_taxon = d.groupby("taxon")["diff"].mean().dropna()
        if len(per_taxon) < min_pairs:
            continue
        nz = per_taxon[per_taxon != 0]
        p_sign = binomtest(int((nz > 0).sum()), len(nz)).pvalue if len(nz) else np.nan
        p_w = wilcoxon(nz).pvalue if len(nz) >= 6 else np.nan
        out.append({"feature": f, "n_taxa": len(per_taxon), "mean_diff": per_taxon.mean(),
                    "frac_positive": float((nz > 0).mean()) if len(nz) else np.nan,
                    "p_sign": p_sign, "p_wilcoxon": p_w})
    res = pd.DataFrame(out)
    if not res.empty:  # Benjamini-Hochberg across all tested features
        p = res["p_sign"].fillna(1.0).values
        order = np.argsort(p)
        q = np.empty_like(p)
        q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        res["q_sign"] = np.minimum(q, 1.0)
    return res.sort_values("p_sign") if not res.empty else res
