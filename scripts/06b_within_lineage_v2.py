"""Single-lineage test with ThermoBase + AciDB added. Subsets: verified literature only vs all tiers; covariates: ph_opt,
and ph_opt + growth temperature. Usage: python3 scripts/06b_within_lineage_v2.py"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
from acidophile.within_lineage import analyse

D = ROOT / "results/within_lineage_v2"
df = (pd.read_csv(D / "lineage_phenotypes.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv"), left_on="ncbi_assembly", right_on="genome_id"))
summ = []
for lineage, g in df.groupby("lineage"):
    for tier, sub in (("verified_only", g[g.tier == "verified"]), ("all_tiers", g)):
        for covs in (("ph_opt",), ("ph_opt", "temp_opt_c")):
            s2 = sub.dropna(subset=list(covs))
            if len(s2) < 6:
                print(f"{lineage}/{tier}/{'+'.join(covs)}: n={len(s2)} too few"); continue
            tab, s = analyse(s2, covariates=covs)
            tab.to_csv(D / f"test_{lineage}_{tier}_{'+'.join(covs)}.csv", index=False)
            summ.append({"lineage": lineage, "subset": tier, **s, "top_feature": tab.iloc[0].feature, "top_rho": tab.iloc[0].rho_with_residual,
                         "top_p_uncorrected": tab.iloc[0].p_uncorrected, "top_p_family_wise": tab.iloc[0].p_family_wise})
S = pd.DataFrame(summ)
S.to_csv(D / "summary.csv", index=False)
cols = ["lineage", "subset", "covariates", "n", "opt_r2", "residual_sd", "detectable_abs_rho_family_wise", "top_feature", "top_rho", "top_p_uncorrected", "top_p_family_wise"]
print(S[cols].round(3).to_string(index=False))
