"""Full-dataset model comparison (M1-M6 with and without pH_opt) under family-blocked and leave-one-lineage-out validation,
on all species, by evidence tier, and by domain.
Usage: python3 scripts/11_full_models.py   (reads results/full/{phenotypes_all,genome_features,hmm_counts}.csv)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
from acidophile.models import compare_models

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/full")      # DATA_DIR=results/broad for the expanded set
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df["tier"] = df.confidence.map({"high": "verified", "medium": "verified", "low": "database_or_unverified"})
df = df[df.ph_opt.notna()].reset_index(drop=True)
print(f"{len(df)} species with genome features and pH_opt | pH_min range {df.ph_min.min():.1f}-{df.ph_min.max():.1f}, sd {df.ph_min.std():.2f} | "
      f"{df.domain.value_counts().to_dict()} | families {df.family.nunique()} | classes {df['class'].nunique()}")
print("pH_min by origin:", df.groupby("origin").ph_min.agg(["count", "mean", "std"]).round(2).to_dict("index"))
subsets = {"all": df, "verified_only": df[df.tier == "verified"], "database_only": df[df.tier != "verified"],
           "bacteria": df[df.domain == "Bacteria"], "archaea": df[df.domain == "Archaea"]}
rows, full = [], []
for sname, sub in subsets.items():
    for scheme, grp in (("family_blocked", "family"), ("leave_one_lineage_out", "class")):
        try:
            res, _ = compare_models(sub, group_col=grp, scheme=scheme)
        except ValueError as e:
            print(f"{sname}/{scheme}: skipped ({e})"); continue
        res.insert(0, "subset", sname); res.insert(1, "scheme", scheme); full.append(res)
        b = res[res.model == "M1_opt_only"].iloc[0]
        r = res[(res.model != "M1_opt_only") & ~res.model.str.contains("noopt")].sort_values("delta_rmse").iloc[0]
        rows.append({"subset": sname, "scheme": scheme, "n": int(b.n), "M1_rmse": round(b.rmse, 3), "M1_r2": round(b.r2, 3),
                     "best_other": r.model, "best_delta_rmse": round(r.delta_rmse, 3), "ci": f"{r.delta_lo:.3f} to {r.delta_hi:.3f}",
                     "models_beating_M1": ", ".join(res[res.beats_baseline == True].model) or "none"})
pd.concat(full).to_csv(D / "model_comparison_all_subsets.csv", index=False)
S = pd.DataFrame(rows); S.to_csv(D / "model_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(S.to_string(index=False))
