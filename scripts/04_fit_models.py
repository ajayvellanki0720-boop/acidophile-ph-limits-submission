"""Usage: python3 scripts/04_fit_models.py phenotypes_primary.csv genome_features.csv system_counts.csv hmm_counts.csv OUTDIR [keyword]
Runs M1-M6 under family-blocked and leave-one-lineage-out CV, overall and by domain / acid class."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.models import compare_by_stratum, compare_models

ph, gf, sc, hmm, out = sys.argv[1:6]
source = sys.argv[6] if len(sys.argv) > 6 else "hmm"
out = Path(out); out.mkdir(parents=True, exist_ok=True)
df = (pd.read_csv(ph).merge(pd.read_csv(gf), on="ncbi_assembly")
      .merge(pd.read_csv(sc), left_on="ncbi_assembly", right_on="genome_id")
      .merge(pd.read_csv(hmm).drop(columns=["n_genes"]), on="genome_id"))
for scheme, grp in [("family_blocked", "family"), ("leave_one_lineage_out", "class")]:
    res, preds = compare_models(df, group_col=grp, scheme=scheme, source=source)
    res.to_csv(out / f"models_{scheme}.csv", index=False); preds.to_csv(out / f"preds_{scheme}.csv", index=False)
    print(f"\n== {scheme} (n={res.n.iloc[0]}) ==\n", res.round(3).to_string(index=False))
    for col in ("domain", "acid_class"):
        for k, v in compare_by_stratum(df, col, group_col=grp, scheme=scheme, source=source).items():
            if isinstance(v, str): print(f"[{col}={k}] skipped: {v}")
            else: v.to_csv(out / f"models_{scheme}_{col}_{k}.csv", index=False)
