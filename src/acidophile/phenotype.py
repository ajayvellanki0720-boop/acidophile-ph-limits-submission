"""Stage 1: phenotype table schema, bound-type handling, and acid chemistry.

pH_min (lowest pH permitting growth) is the primary outcome. A value is only an
exact limit if growth was shown to fail below it. "Grew at pH 1.0, lower not
tested" is a censored bound and is kept out of the primary set.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

BOUND_TYPES = {"measured", "lowest_tested", "highest_tested", "extrapolated", "unknown"}
CONFIDENCE = {"high", "medium", "low"}
ACID_IDENTITIES = {"H2SO4", "HCl", "organic", "mixed", "buffered_unspecified", "unspecified"}

REQUIRED_COLUMNS = [
    "strain_id", "ncbi_assembly", "domain", "phylum", "class", "family", "genus", "species",
    "ph_min", "ph_min_bound", "ph_opt", "ph_max", "ph_max_bound",
    "acid_identity", "organic_acid", "organic_acid_mM", "ph_measured_at_c",
    "temp_opt_c", "temp_min_c", "temp_max_c", "oxygen", "medium", "metabolism",
    "metal_note", "source", "confidence",
]


def proton_activity(ph):
    """a_H+ = 10^-pH (molar activity). Growth limits are compared on pH, not a_H."""
    return np.power(10.0, -np.asarray(ph, dtype=float))


def undissociated_fraction(ph, pka):
    """Fraction of a monoprotic weak acid in the membrane-permeant HA form.

    Equal pH with different acids is not equal stress: acetic acid (pKa 4.76)
    is ~64% HA at pH 4.5 and uncouples the proton gradient, while sulfate is
    essentially never permeant. Used to weight organic-acid media.
    """
    return 1.0 / (1.0 + np.power(10.0, np.asarray(ph, dtype=float) - pka))


def classify_acid_class(ph_opt, ph_max, ph_min, opt_cut=5.0, extreme_cut=3.0, obligate_max=6.5):
    """extreme / moderate acidophile, acid_tolerant, or non_acidophile."""
    if pd.isna(ph_opt):
        return "unknown"
    if ph_opt <= extreme_cut:
        return "extreme_acidophile"
    if ph_opt <= opt_cut:
        return "moderate_acidophile"
    if not pd.isna(ph_min) and ph_min <= 4.5:
        return "acid_tolerant"
    return "non_acidophile"


def validate(df: pd.DataFrame):
    """Return (annotated_df, issues_df). Never silently repairs values."""
    issues = []

    def flag(i, field, problem):
        issues.append({"strain_id": df.at[i, "strain_id"] if "strain_id" in df else i,
                       "field": field, "problem": problem})

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"missing required columns: {missing}")
    out = df.copy()
    for c in ["ph_min", "ph_opt", "ph_max", "temp_opt_c", "temp_min_c", "temp_max_c",
              "organic_acid_mM", "ph_measured_at_c"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    if out["strain_id"].duplicated().any():
        for i in out.index[out["strain_id"].duplicated(keep=False)]:
            flag(i, "strain_id", "duplicate strain_id")
    for i, r in out.iterrows():
        for c in ("ph_min", "ph_opt", "ph_max"):
            if not pd.isna(r[c]) and not (-1.0 <= r[c] <= 14.0):
                flag(i, c, f"{c}={r[c]} outside plausible -1..14")
        vals = [r["ph_min"], r["ph_opt"], r["ph_max"]]
        present = [v for v in vals if not pd.isna(v)]
        if present != sorted(present):
            flag(i, "ph_min/ph_opt/ph_max", "not ordered min <= opt <= max")
        for c in ("ph_min_bound", "ph_max_bound"):
            if r[c] not in BOUND_TYPES:
                flag(i, c, f"bound type {r[c]!r} not in {sorted(BOUND_TYPES)}")
        if r["confidence"] not in CONFIDENCE:
            flag(i, "confidence", f"{r['confidence']!r} not high/medium/low")
        if r["acid_identity"] not in ACID_IDENTITIES:
            flag(i, "acid_identity", f"{r['acid_identity']!r} not recognised")
        if r["acid_identity"] in ("organic", "mixed") and pd.isna(r["organic_acid_mM"]):
            flag(i, "organic_acid_mM", "organic acid present but concentration missing")
        if pd.isna(r["ph_measured_at_c"]):
            flag(i, "ph_measured_at_c", "temperature of pH measurement missing (matters for thermoacidophiles)")
        if pd.isna(r["ncbi_assembly"]) or not str(r["ncbi_assembly"]).startswith(("GCA_", "GCF_")):
            flag(i, "ncbi_assembly", "no exact matched NCBI assembly")

    out["acid_class"] = [classify_acid_class(r.ph_opt, r.ph_max, r.ph_min) for r in out.itertuples()]
    out["acid_reach"] = out["ph_opt"] - out["ph_min"]      # pH units tolerated below the optimum
    out["ph_breadth"] = out["ph_max"] - out["ph_min"]      # specialist (narrow) vs generalist (wide)
    out["aH_min"] = proton_activity(out["ph_min"])
    out["ph_min_is_exact"] = out["ph_min_bound"].eq("measured") & out["ph_min"].notna()
    return out, pd.DataFrame(issues, columns=["strain_id", "field", "problem"])


def primary_set(df: pd.DataFrame, min_confidence=("high", "medium")) -> pd.DataFrame:
    """Rows eligible for the primary pH_min model: exact bound, matched genome, ph_opt known."""
    m = (df["ph_min_is_exact"] & df["ph_opt"].notna() & df["confidence"].isin(min_confidence)
         & df["ncbi_assembly"].astype(str).str.startswith(("GCA_", "GCF_")))
    return df[m].copy()


def censored_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Rows where ph_min is only 'lowest tested': true pH_min <= reported value."""
    return df[df["ph_min_bound"].eq("lowest_tested") & df["ph_min"].notna()].copy()
