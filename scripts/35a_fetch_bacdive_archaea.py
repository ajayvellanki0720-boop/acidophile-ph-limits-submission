"""Fetch raw BacDive records for EVERY archaeal genus listed by NCBI Taxonomy -> data/raw/bacdive_archaea/<genus>.json (public API, polite delays, cached; raw files never edited).
Genus list: NCBI Datasets filtered subtree of Archaea (taxid 2157) at rank GENUS, saved to data/curated/archaeal_genera.txt. Usage: python3 scripts/35a_fetch_bacdive_archaea.py"""
import importlib.util
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bd", ROOT / "scripts/00_fetch_bacdive.py"); bd = importlib.util.module_from_spec(spec); spec.loader.exec_module(bd)
OUT = ROOT / "data/raw/bacdive_archaea"; OUT.mkdir(parents=True, exist_ok=True)
GEN = ROOT / "data/curated/archaeal_genera.txt"
if not GEN.exists():
    d = json.load(urllib.request.urlopen("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/2157/filtered_subtree?specified_limit=false&rank_limits=GENUS", timeout=90))
    GEN.write_text("\n".join(sorted(v["scientific_name"] for v in d["edges"].values() if v.get("rank") == "GENUS" and " " not in v["scientific_name"])) + "\n")
genera = GEN.read_text().split()
print(len(genera), "archaeal genera", flush=True)
n = 0
for g in genera:
    f = OUT / f"{g}.json"
    if f.exists():
        continue
    try:
        recs = bd.fetch_genus(g)
    except Exception as e:
        print(g, "no BacDive records or error:", str(e)[:60], flush=True); f.write_text("{}"); continue
    f.write_text(json.dumps(recs)); n += 1
    print(g, len(recs), flush=True)
print("DONE", n, "genera with records", flush=True)
