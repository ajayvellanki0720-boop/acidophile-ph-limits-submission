"""Usage: python3 scripts/06_within_lineage.py   (reads results/within_lineage/*, writes results/within_lineage/test_*.csv)"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
from acidophile.within_lineage import analyse

D = ROOT / "results/within_lineage"
df = (pd.read_csv(D / "lineage_phenotypes.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv"), left_on="ncbi_assembly", right_on="genome_id"))
for lineage, g in df.groupby("lineage"):
    for tier, sub in (("all", g), ("medium_or_high_only", g[g.confidence.isin(["medium", "high"])])):
        if len(sub) < 6:
            print(f"{lineage}/{tier}: n={len(sub)} too few"); continue
        tab, s = analyse(sub)
        tab.to_csv(D / f"test_{lineage}_{tier}.csv", index=False)
        print(f"\n=== {lineage} / {tier} ===")
        print({k: round(v, 3) if isinstance(v, float) else v for k, v in s.items()})
        print(tab.head(6).round(3).to_string(index=False))
