"""Resolve species -> NCBI assembly (versioned accession) via the public NCBI Datasets API.
Saves every candidate assembly to data/raw/ncbi_assemblies.json (audit trail), then picks one per species:
  1) assembly from type material, 2) RefSeq reference/representative, 3) complete genome, 4) best contig N50.
Strain-level matching to the phenotype strain is a separate, explicit check (column `type_material`).
Usage: python3 scripts/00e_resolve_assemblies.py data/curated/taxa_to_resolve.txt data/curated/assembly_choice.csv"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report?page_size=100&filters.assembly_source=all"
LEVEL = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}


def reports(taxon):
    try:
        d = json.load(urllib.request.urlopen(API.format(urllib.parse.quote(taxon)), timeout=60))
    except Exception:
        return []
    return d.get("reports", [])


def summarise(r):
    ai, org = r.get("assembly_info", {}), r.get("organism", {})
    st = r.get("assembly_stats", {})
    return {"accession": r["accession"], "organism": org.get("organism_name"),
            "tax_id": org.get("tax_id"), "strain": org.get("infraspecific_names", {}).get("strain"),
            "level": ai.get("assembly_level"), "refseq_category": ai.get("refseq_category"),
            "type_material": (r.get("type_material") or {}).get("type_display_text"),
            "contig_n50": st.get("contig_n50", 0), "total_length": int(st.get("total_sequence_length", 0)),
            "checkm_completeness": (r.get("checkm_info") or {}).get("completeness"),
            "source_database": r.get("source_database")}


EPITHET_ALIAS = {"organivorus": "organivorans"}  # orthographic corrections only


def epithet_match(taxon, organism):
    """Same species epithet (genus may differ: reclassification). A different epithet is a synonym or another
    species, whose strain is NOT the strain the phenotype describes, so it is rejected."""
    want = EPITHET_ALIAS.get(taxon.split()[1], taxon.split()[1])
    return isinstance(organism, str) and len(organism.split()) > 1 and organism.split()[1] == want


def rank(c):
    return (c["type_material"] is not None, c["refseq_category"] in ("reference genome", "representative genome"),
            LEVEL.get(c["level"], 0), c["contig_n50"] or 0)


if __name__ == "__main__":
    taxa = [t.strip() for t in Path(sys.argv[1]).read_text().splitlines() if t.strip()]
    raw, rows = {}, []
    for t in taxa:
        binom = " ".join(t.split()[:2])
        allc = [summarise(r) for r in reports(binom)]
        raw[t] = allc
        time.sleep(0.25)
        cands = [c for c in allc if epithet_match(t, c["organism"])]
        if not cands:
            rows.append({"taxon": t, "n_assemblies": 0, "n_rejected_other_species": len(allc)}); continue
        best = max(cands, key=rank)
        rows.append({"taxon": t, "n_assemblies": len(cands), **best})
    rawf = ROOT / "data/raw/ncbi_assemblies.json"
    old = json.loads(rawf.read_text()) if rawf.exists() else {}
    old.update(raw)                      # merge: a partial run must never discard earlier candidates
    rawf.write_text(json.dumps(old))
    pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
    df = pd.DataFrame(rows)
    print(f"{len(taxa)} taxa | {(df.n_assemblies > 0).sum()} with assemblies | {df.type_material.notna().sum()} with type-material assembly")
