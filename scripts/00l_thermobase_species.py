"""Collapse ThermoBase strain-level rows to one record per species (median), keeping acid-relevant species
(ph_min <= 5.5 or ph_opt <= 6.5). Species whose rows disagree by >1 pH unit in ph_min are flagged `conflict`.
Output: data/curated/thermobase_species.csv. Bounds stay `unknown`: ThermoBase does not say whether a minimum is a limit."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
b = pd.read_csv(ROOT / "data/curated/thermobase_candidates.csv")
b = b[b.species.notna() & b.ph_min.notna() & ((b.ph_min <= 5.5) | (b.ph_opt <= 6.5))]
agg = b.groupby("species").agg(n_rows=("Name", "size"), taxid=("taxid", "first"), domain=("Domain", "first"), phylum=("Phylum", "first"),
                               order=("Order", "first"), family=("Family", "first"),
                               ph_min=("ph_min", "median"), ph_opt=("ph_opt", "median"), ph_max=("ph_max", "median"),
                               ph_min_spread=("ph_min", lambda s: s.max() - s.min()), temp_opt_c=("temp_opt_c", "median"),
                               sources=("Source", lambda s: " | ".join(map(str, s.dropna().unique()[:3])))).reset_index()
agg["conflict"] = agg.ph_min_spread > 1.0
agg.to_csv(ROOT / "data/curated/thermobase_species.csv", index=False)
print(len(agg), "species |", int(agg.conflict.sum()), "with conflicting rows |", int((agg.n_rows > 1).sum()), "multi-row")
cur = pd.read_csv(ROOT / "data/curated/taxa_to_resolve.txt", header=None)[0].tolist()
allt = sorted(set(cur) | set(agg.species))
Path(ROOT / "data/curated/taxa_to_resolve.txt").write_text("\n".join(allt))
print(len(allt), "taxa to resolve")
