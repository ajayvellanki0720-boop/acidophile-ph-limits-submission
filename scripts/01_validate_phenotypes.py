"""Validate a curated phenotype CSV (see data/templates/phenotype_template.csv).
Usage: python3 scripts/01_validate_phenotypes.py data/curated/phenotypes.csv
Writes data/curated/phenotypes_validated.csv, phenotype_issues.csv, phenotypes_primary.csv, phenotypes_censored.csv."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.phenotype import censored_rows, primary_set, validate

src = Path(sys.argv[1])
out_dir = src.parent
df, issues = validate(pd.read_csv(src))
df.to_csv(out_dir / "phenotypes_validated.csv", index=False)
issues.to_csv(out_dir / "phenotype_issues.csv", index=False)
primary_set(df).to_csv(out_dir / "phenotypes_primary.csv", index=False)
censored_rows(df).to_csv(out_dir / "phenotypes_censored.csv", index=False)
print(f"{len(df)} strains | {len(issues)} issues | primary {len(primary_set(df))} | censored {len(censored_rows(df))}")
print(df["acid_class"].value_counts().to_string())
