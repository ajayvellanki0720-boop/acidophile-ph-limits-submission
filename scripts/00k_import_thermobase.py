"""Import pH fields from ThermoBase (read in place, not copied) and Topt from TEMPURA (temperature only, no pH).
Writes data/curated/thermobase_candidates.csv and data/curated/tempura_topt.csv.
ThermoBase does not define whether Min/Max pH are measured limits or tested edges, so every bound is `unknown`
(confidence low) until checked against the record's own Source link.
Usage: python3 scripts/00k_import_thermobase.py /path/to/ThermoBase.xlsx /path/to/200617_TEMPURA.csv"""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
tb = pd.read_excel(sys.argv[1])
tb.columns = [c.strip() for c in tb.columns]
rows = []
for r in tb.itertuples(index=False):
    d = r._asdict() if hasattr(r, "_asdict") else {}
tb = tb.rename(columns={"Min. pH": "ph_min", "Max. pH": "ph_max", "Avg. Optimum pH": "ph_opt", "Avg. Optimum Temp.": "temp_opt_c",
                        "Min. Temp. (°C)": "temp_min_c", "Max. Temp. (°C)": "temp_max_c", "Taxonomic ID": "taxid"})
tb["species"] = tb["Name"].str.replace(r'["“”]', "", regex=True).str.extract(r"^(Candidatus\s+)?([A-Z][a-z]+ [a-z]+)")[1]
tb["strain_text"] = tb["Name"].str.replace(r'^(Candidatus\s+)?[A-Z]\w+ \w+\s*', "", regex=True)
def _clean(v):
    """ThermoBase has transcription typos such as '3 .5' and '7.7.'; repaired and flagged, never silently."""
    if pd.isna(v) or isinstance(v, (int, float)):
        return v, False
    s = str(v).replace(" ", "").rstrip(".")
    try:
        return float(s), True
    except ValueError:
        return float("nan"), True


for col in ("ph_min", "ph_opt", "ph_max"):
    tb[col], fixed = zip(*[_clean(v) for v in tb[col]])
    tb[col + "_typo_repaired"] = fixed
out = tb[tb.ph_min.notna() | tb.ph_opt.notna()][["Name", "species", "strain_text", "taxid", "Domain", "Phylum", "Class", "Order", "Family",
                                                   "ph_min", "ph_opt", "ph_max", "ph_min_typo_repaired", "ph_opt_typo_repaired", "ph_max_typo_repaired", "temp_min_c", "temp_opt_c", "temp_max_c", "Oxygen Requirement", "Source"]]
out.to_csv(ROOT / "data/curated/thermobase_candidates.csv", index=False)
t = pd.read_csv(sys.argv[2])
t = t.rename(columns={"genus_and_species": "species", "Topt_ave": "tempura_topt_c", "Tmin": "tempura_tmin_c", "Tmax": "tempura_tmax_c"})
t["species"] = t.species.str.replace('"', "", regex=False).str.split().str[:2].str.join(" ")
t[["species", "taxonomy_id", "strain", "tempura_tmin_c", "tempura_topt_c", "tempura_tmax_c"]].to_csv(ROOT / "data/curated/tempura_topt.csv", index=False)
print(f"ThermoBase: {len(out)} rows with pH ({out.ph_min.notna().sum()} min, {out.ph_opt.notna().sum()} opt, {out.ph_max.notna().sum()} max); TEMPURA: {len(t)} rows")
