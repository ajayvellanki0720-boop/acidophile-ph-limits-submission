"""Candidate NEW archaeal species for the enlarged archaeal sample -> data/curated/archaea_ext_candidates.csv and taxa list data/curated/taxa_archaea_ext.txt.
Sources, in precedence order: (1) parser-extracted species-description abstracts (broad_candidates.csv, rows flagged clean or usable_noopt) for taxa whose genus is archaeal
and whose species are not already in results/broad2/phenotypes_all.csv; (2) BacDive archaeal strains (35b) with a growth range >= 1 pH unit (median over strains).
Usage: python3 scripts/35c_archaea_candidates.py"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ph = pd.read_csv(ROOT / "results/broad2/phenotypes_all.csv")
have = set(ph[ph.domain == "Archaea"].species.astype(str)) | set(ph.species.astype(str))
arch_genera = set(ph[ph.domain == "Archaea"].genus) | set((ROOT / "data/curated/archaeal_genera.txt").read_text().split())
b = pd.read_csv(ROOT / "data/curated/broad_candidates.csv")
b = b[(b.clean | b.usable_noopt) & b.taxon.astype(str).str.split().str[0].isin(arch_genera)]
b = b[~b.taxon.isin(have)].copy()
b["species"] = b.taxon.str.split().str[:2].str.join(" ")
lit = (b.sort_values("clean", ascending=False).groupby("species").first().reset_index()[["species", "ph_min", "ph_opt", "ph_max", "pmid", "doi"]]
       .assign(origin="literature_auto", source=lambda d: "PMID " + d.pmid.astype("Int64").astype(str)))
c = pd.read_csv(ROOT / "data/curated/bacdive_archaea_candidates.csv")
c = c[c.ph_min.notna() & c.ph_max.notna() & (c.ph_max - c.ph_min >= 1.0) & ~c.species.isin(have)]
bd = (c.groupby("species").agg(ph_min=("ph_min", "median"), ph_opt=("ph_opt", "median"), ph_max=("ph_max", "median"), source=("source", "first"), temp_opt_c=("temp_opt_c", "median"),
                              phylum=("phylum", "first"), cls=("class", "first"), family=("family", "first"), genus=("genus", "first")).reset_index().assign(origin="bacdive"))
bd = bd[~bd.species.isin(lit.species)]
out = pd.concat([lit, bd], ignore_index=True)
out.to_csv(ROOT / "data/curated/archaea_ext_candidates.csv", index=False)
(ROOT / "data/curated/taxa_archaea_ext.txt").write_text("\n".join(out.species) + "\n")
print(len(lit), "from abstracts +", len(bd), "from BacDive =", len(out), "new archaeal species")
