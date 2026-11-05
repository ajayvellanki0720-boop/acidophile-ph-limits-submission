"""Assemble data/curated/phenotypes.csv (template schema) from literature_curated.csv, bacdive_candidates.csv,
assembly_choice.csv and NCBI taxonomy. Rules:
  - literature value wins over BacDive for the same species; disagreements >0.5 pH are logged
  - confidence `medium` only for a literature `measured` bound whose assembly is from type material and whose
    strain assignment is not an assumption; everything else is `low` (kept, but outside the primary set)
  - species with no epithet-matched assembly are dropped here (nothing to model), and counted
Usage: python3 scripts/00g_build_phenotype_table.py"""
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from acidophile.phenotype import REQUIRED_COLUMNS

cur = pd.read_csv(ROOT / "data/curated/literature_curated.csv", dtype={"pmid": str})
bd = pd.read_csv(ROOT / "data/curated/bacdive_candidates.csv")
asm = pd.read_csv(ROOT / "data/curated/assembly_choice.csv").set_index("taxon")

rows = []
for r in cur.itertuples():
    rows.append(dict(taxon=r.taxon, ph_min=r.ph_min, ph_opt=r.ph_opt, ph_max=r.ph_max,
                     bound="lowest_tested" if r.bound == "censored" else "measured", src=f"PMID:{r.pmid}", origin="literature",
                     assumption=bool(r.assumption), temp_opt_c=r.temp_opt_c, note=r.note, min_verified=bool(r.min_verified)))
lit_taxa = {x["taxon"] for x in rows}
disc = []
# precedence: literature > ThermoBase > BacDive. ThermoBase bounds are unverified (it does not define them).
tbs = pd.read_csv(ROOT / "data/curated/thermobase_species.csv")
tb_taxa = set()
for r in tbs.itertuples():
    if r.conflict or r.species in lit_taxa:
        if r.species in lit_taxa:
            l = next(x for x in rows if x["taxon"] == r.species)
            if pd.notna(l["ph_min"]) and abs(r.ph_min - l["ph_min"]) > 0.5:
                disc.append({"taxon": r.species, "lit_ph_min": l["ph_min"], "other_ph_min": r.ph_min, "lit_src": l["src"], "other": "ThermoBase"})
        continue
    tb_taxa.add(r.species)
    rows.append(dict(taxon=r.species, ph_min=r.ph_min, ph_opt=r.ph_opt, ph_max=r.ph_max, bound="unknown", src="ThermoBase:" + (re.search(r"\d+", str(r.taxid)) or ["NA"])[0],
                     origin="thermobase", assumption=False, temp_opt_c=r.temp_opt_c, note=f"ThermoBase median of {r.n_rows} row(s); bound type unverified", min_verified=False))
# AciDB (https://acidb.cl): genome-linked acidophile table. Only rows with a real growth range and a non-derived optimum
# (isolated, ph_min < ph_max, optimum not the range midpoint, assembly quality ok) are used; one genome per species.
ac = pd.read_csv(ROOT / "data/curated/acidb_organisms.csv")
ac = ac[ac.usable].copy()
ac["q"] = (ac.state.eq("complete")).astype(int) * 1000 + ac.gen_completeness - 5 * ac.gen_contamination
ac = ac.sort_values("q", ascending=False)
epi = lambda s: s.split()[1] if len(s.split()) > 1 else s
have_epi = {epi(x) for x in lit_taxa | tb_taxa}
acidb_override = {}
for sp, g in ac.groupby("species", sort=False):
    if sp in lit_taxa | tb_taxa or epi(sp) in have_epi:
        continue
    best = g.iloc[0]
    rows.append(dict(taxon=sp, ph_min=g.ph_min.median(), ph_opt=g.ph_associated.median(), ph_max=g.ph_max.median(), bound="unknown",
                     src=f"AciDB:{int(best.id_organism)}", origin="acidb", assumption=False, temp_opt_c=g.temp_associated.median(),
                     note=f"AciDB median of {len(g)} genome row(s); bound type unverified", min_verified=False))
    acidb_override[sp] = (best.access_id, best.strain_names, best.state)
# Parser-extracted species descriptions (scripts 22-23; audited sample: 40 of 40 rows consistent with the abstract). Added last among literature
# sources and only for species not already present from any other source.
bf = ROOT / "data/curated/broad_with_genomes.csv"
if bf.exists():
    bw = pd.read_csv(bf, dtype={"pmid": str})
    have_now = {x["taxon"] for x in rows}
    for r in bw.itertuples():
        if r.taxon in have_now or epi(r.taxon) in {epi(x) for x in have_now if len(x.split()) > 1 and x.split()[0] == r.taxon.split()[0]}:
            continue
        rows.append(dict(taxon=r.taxon, ph_min=r.ph_min, ph_opt=r.ph_opt, ph_max=r.ph_max, bound="measured", src=f"PMID:{r.pmid}", origin="literature_auto",
                         assumption=False, temp_opt_c=np.nan, note="parser-extracted from abstract (single-taxon title, no multi-strain or hedged wording)", min_verified=False))
lit_taxa_all = lit_taxa | tb_taxa | set(acidb_override)
for r in bd.itertuples():
    sp = r.species
    if sp in lit_taxa_all:
        l = next(x for x in rows if x["taxon"] == sp)
        if pd.notna(r.ph_min) and pd.notna(l["ph_min"]) and abs(r.ph_min - l["ph_min"]) > 0.5:
            disc.append({"taxon": sp, "lit_ph_min": l["ph_min"], "bacdive_ph_min": r.ph_min, "lit_src": l["src"], "bacdive": r.source})
        continue
    if pd.isna(r.ph_min):
        continue
    rows.append(dict(taxon=sp, ph_min=r.ph_min, ph_opt=r.ph_opt, ph_max=r.ph_max, bound="unknown", src=r.source,
                     origin="bacdive", assumption=False, min_verified=False, temp_opt_c=r.temp_opt_c, note="BacDive range, bound type unverified"))
df = pd.DataFrame(rows)
df = df[~df.taxon.str.split().str[1].isin(["sp.", "sp", "bacterium", "archaeon"])]      # unnamed species cannot be matched to a type strain
df["n_in"] = len(df)
df = df.join(asm, on="taxon", rsuffix="_asm")
for sp, (acc, strains, state) in acidb_override.items():     # AciDB gives the exact genome, so skip name-based resolution
    i = df.index[(df.taxon == sp) & (df.origin == "acidb")]
    df.loc[i, "accession"] = acc
    df.loc[i, "strain"] = strains
    df.loc[i, "level"] = "Complete Genome" if state == "complete" else "Scaffold"
dropped = df[df.accession.isna()]
df = df[df.accession.notna()].copy()

# taxonomy lineage from NCBI
tax = {}
ids = sorted({str(int(x)) for x in df.tax_id.dropna()})  # AciDB rows have no tax_id yet
for i in range(0, len(ids), 50):
    d = json.load(urllib.request.urlopen("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(ids[i:i + 50]) + "/dataset_report", timeout=60))
    for rep in d.get("reports", []):
        if "taxonomy" in rep:
            c = rep["taxonomy"]["classification"]
            tax[str(rep["taxonomy"]["tax_id"])] = {k: c.get(k, {}).get("name") for k in ("domain", "kingdom", "phylum", "class", "order", "family", "genus")}
    time.sleep(0.3)
# rows without a tax_id (AciDB): look the lineage up by species name
import urllib.parse
by_name = {}
for sp in sorted(set(df.loc[df.tax_id.isna(), "taxon"])):
    try:
        d = json.load(urllib.request.urlopen("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + urllib.parse.quote(sp) + "/dataset_report", timeout=60))
        rep = d["reports"][0]["taxonomy"]
        c = rep["classification"]
        by_name[sp] = {k: c.get(k, {}).get("name") for k in ("domain", "kingdom", "phylum", "class", "order", "family", "genus")}
        df.loc[(df.taxon == sp) & df.tax_id.isna(), "tax_id"] = rep["tax_id"]
    except Exception:
        by_name[sp] = {}
    time.sleep(0.3)
L = df.apply(lambda r: by_name.get(r.taxon) or (tax.get(str(int(r.tax_id)), {}) if pd.notna(r.tax_id) else {}), axis=1)
for k in ("domain", "phylum", "class", "order", "family", "genus"):
    df[k] = L.map(lambda d, k=k: d.get(k))
df["domain"] = df["domain"].fillna(L.map(lambda d: d.get("kingdom")))
df["family_missing"] = df.family.isna()
df["family"] = df.family.fillna(df.genus)     # block at genus where NCBI has no family; flagged
df["class"] = df["class"].fillna(df.phylum)

# Assemblies NCBI has not flagged as type material but whose strain designation equals the protologue's type strain (checked by hand).
TYPE_STRAIN_CHECKED = {"Alicyclobacillus curvatus": "ALEF1"}      # protologue: type strain ALEF1T = CGMCC 1.17055T (PMC9558579)
strain_ok = pd.Series([TYPE_STRAIN_CHECKED.get(tx) == str(st) for tx, st in zip(df.taxon, df.strain)], index=df.index)
ok = (df.bound.eq("measured") & df.origin.isin(["literature", "literature_auto"]) & (df.type_material.notna() | strain_ok) & ~df.assumption)
df["confidence"] = ok.map({True: "medium", False: "low"})
df.loc[ok & df.min_verified, "confidence"] = "high"   # tested range extended beyond the observed minimum, or explicit no-growth statement
out = pd.DataFrame({
    "strain_id": df.taxon.str.replace(" ", "_") + "|" + df.accession, "ncbi_assembly": df.accession,
    "domain": df.domain, "phylum": df.phylum, "class": df["class"], "family": df.family, "genus": df.genus, "species": df.taxon,
    "ph_min": df.ph_min, "ph_min_bound": df.bound, "ph_opt": df.ph_opt, "ph_max": df.ph_max,
    "ph_max_bound": df.bound.map({"measured": "measured", "lowest_tested": "highest_tested", "unknown": "unknown"}),
    "acid_identity": "unspecified", "organic_acid": None, "organic_acid_mM": None, "ph_measured_at_c": None,
    "temp_opt_c": df.temp_opt_c, "temp_min_c": None, "temp_max_c": None, "oxygen": None, "medium": None,
    "metabolism": None, "metal_note": None, "source": df.src, "confidence": df.confidence,
})
tem = pd.read_csv(ROOT / "data/curated/tempura_topt.csv").groupby("species").tempura_topt_c.median()
out["temp_source"] = None
have_t = out.temp_opt_c.notna()
out.loc[have_t, "temp_source"] = "abstract_or_thermobase"
fill = out.temp_opt_c.isna() & out.species.isin(tem.index)
out.loc[fill, "temp_opt_c"] = out.loc[fill, "species"].map(tem)
out.loc[fill, "temp_source"] = "TEMPURA"
out["origin"] = df.origin.values
out["assembly_strain"] = df.strain.values
out["assembly_type_material"] = df.type_material.notna().values
out["assembly_level"] = df.level.values
out["family_missing"] = df.family_missing.values
out["note"] = df.note.values
for c in REQUIRED_COLUMNS:
    assert c in out
# One row per genome: the same organism can enter under two names (e.g. Sulfolobus/Saccharolobus solfataricus).
# Same-source duplicates are merged by median; across sources the higher-precedence source wins.
rank = {"literature": 0, "literature_auto": 0.5, "thermobase": 1, "acidb": 2, "bacdive": 3}
out["_rank"] = out.origin.map(rank)
num = ["ph_min", "ph_opt", "ph_max", "temp_opt_c"]
merged = []
for acc, g in out.sort_values("_rank").groupby("ncbi_assembly", sort=False):
    best = g[g._rank == g._rank.min()]
    row = best.iloc[0].copy()
    if len(best) > 1:
        row[num] = best[num].median()
        row["note"] = f"{row['note']}; merged {len(best)} same-source rows"
    row["n_names_merged"] = len(g)
    merged.append(row)
out = pd.DataFrame(merged).drop(columns="_rank").reset_index(drop=True)
out = out[out.domain.isin(["Bacteria", "Archaea"])]      # AciDB lists a few eukaryotes; out of scope
out = out.drop_duplicates(["species", "source"]).reset_index(drop=True)
out.to_csv(ROOT / "data/curated/phenotypes.csv", index=False)
pd.DataFrame(disc).to_csv(ROOT / "data/curated/source_discrepancies.csv", index=False)
print(f"{len(out)} rows with matched assembly | dropped (no species-matched assembly): {len(dropped)} | "
      f"discrepancies BacDive vs literature: {len(disc)}")
print(out.confidence.value_counts().to_dict(), out.domain.value_counts().to_dict(), "| families:", out.family.nunique())
