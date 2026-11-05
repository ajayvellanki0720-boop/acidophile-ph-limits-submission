"""Enlarged archaeal table -> results/archaea_ext/phenotypes_all.csv: the existing archaeal rows of results/broad2/phenotypes_all.csv plus the new species (35c) that resolved to an assembly
(35c + 00e). Taxonomy for new rows from NCBI Taxonomy by tax_id; assemblies need CheckM completeness >= 90 % where reported. New rows are flagged `ext` = True.
Optional arguments: candidates.csv assembly_choice.csv base_phenotypes.csv outdir (defaults: the first extension). Usage: python3 scripts/35d_build_archaea_ext_table.py [cand assemblies base outdir]"""
import json, sys, time, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARG = sys.argv[1:]
OUT = ROOT / (ARG[3] if len(ARG) > 3 else "results/archaea_ext"); OUT.mkdir(parents=True, exist_ok=True)
old = pd.read_csv(ROOT / (ARG[2] if len(ARG) > 2 else "results/broad2/phenotypes_all.csv")); old = old[old.domain == "Archaea"]
old["ext"] = old["ext"] if "ext" in old else False
old["completeness_lt90"] = old["completeness_lt90"] if "completeness_lt90" in old else False
c = pd.read_csv(ROOT / (ARG[0] if ARG else "data/curated/archaea_ext_candidates.csv"))
a = pd.read_csv(ROOT / (ARG[1] if len(ARG) > 1 else "data/curated/assembly_choice_archaea_ext.csv")); a = a[a.n_assemblies > 0]
MINC = float(__import__("os").environ.get("MIN_COMPLETENESS", 90))  # assemblies below 90 % are kept only when MIN_COMPLETENESS is lowered, flagged completeness_lt90
a = a[a.checkm_completeness.isna() | (a.checkm_completeness >= MINC)]
m = c.merge(a, left_on="species", right_on="taxon")
m = m[~m.accession.isin(set(old.ncbi_assembly))]


def lineage(tid):
    r = json.load(urllib.request.urlopen(f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{int(tid)}/dataset_report", timeout=60))["reports"][0]["taxonomy"]
    cl = r.get("classification", {}); g = lambda k: cl.get(k, {}).get("name")
    return g("phylum"), g("class"), g("family"), g("genus")


rows = []
for _, r in m.iterrows():
    try: phy, cls, fam, gen = lineage(r.tax_id)
    except Exception: phy = cls = fam = gen = None
    time.sleep(0.2)
    conf = "medium" if (r.origin in ("literature_auto", "literature") and isinstance(r.type_material, str)) else "low"
    rows.append({"strain_id": f"{r.species.replace(' ', '_')}|{r.accession}", "ncbi_assembly": r.accession, "domain": "Archaea", "phylum": phy or r.get("phylum"), "class": cls or r.get("cls"),
                 "family": fam or r.get("family"), "genus": gen or r.species.split()[0], "species": r.species, "ph_min": r.ph_min, "ph_min_bound": "unknown", "ph_opt": r.ph_opt, "ph_max": r.ph_max,
                 "ph_max_bound": "unknown", "acid_identity": "unspecified", "temp_opt_c": r.get("temp_opt_c"), "source": r.source, "confidence": conf, "origin": r.origin,
                 "assembly_strain": r.strain, "pmid": r.get("pmid"), "assembly_type_material": isinstance(r.type_material, str), "assembly_level": r.level, "note": (r.get("note") if isinstance(r.get("note"), str) else "enlarged-archaea extension") + "; bound type unverified", "ext": True, "completeness_lt90": bool(r.checkm_completeness < 90) if r.checkm_completeness == r.checkm_completeness else False})
new = pd.DataFrame(rows)
allp = pd.concat([old, new], ignore_index=True)
allp.to_csv(OUT / "phenotypes_all.csv", index=False)
print(len(old), "existing +", len(new), "new =", len(allp), "| new classes:", new["class"].value_counts().to_dict(), "| new halophiles:", int((new["class"] == "Halobacteria").sum()))
new[["ncbi_assembly"]].to_csv(OUT / "new_accessions.csv", index=False)
