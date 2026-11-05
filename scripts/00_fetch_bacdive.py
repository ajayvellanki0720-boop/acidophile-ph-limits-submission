"""Fetch raw BacDive records for the candidate acidophile genera into data/raw/bacdive/<genus>.json.
Public API, no login at time of writing; polite delay between calls. Raw files are never edited.
Usage: python3 scripts/00_fetch_bacdive.py [Genus ...]"""
import json
import sys
import time
import urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "bacdive"
GENERA = [
    # bacterial iron/sulfur oxidisers, acidophilic heterotrophs
    "Acidithiobacillus", "Leptospirillum", "Acidiphilium", "Acidocella", "Acidisphaera", "Acidobacterium",
    "Sulfobacillus", "Alicyclobacillus", "Acidimicrobium", "Ferrimicrobium", "Acidiferrobacter", "Thiomonas",
    "Acidibacillus", "Acidihalobacter", "Ferrithrix", "Acidithrix",
    # thermoacidophilic archaea
    "Sulfolobus", "Saccharolobus", "Acidianus", "Metallosphaera", "Sulfurisphaera", "Picrophilus",
    "Thermoplasma", "Ferroplasma", "Acidilobus", "Stygiolobus",
    # acid-tolerant contrasts
    "Acetobacter", "Gluconobacter", "Acidomonas", "Helicobacter", "Streptococcus",
]


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "acidophile-ph-limits research"}), timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def fetch_genus(genus):
    ids = []
    page = f"https://api.bacdive.dsmz.de/taxon/{genus}"
    while page:
        d = get(page)
        r = d["results"]
        ids += list(r.keys() if isinstance(r, dict) else r)
        page = d.get("next")
        time.sleep(0.5)
    recs = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        recs.update(get("https://api.bacdive.dsmz.de/fetch/" + ";".join(map(str, chunk)))["results"])
        time.sleep(0.5)
    return recs


if __name__ == "__main__":
    RAW.mkdir(parents=True, exist_ok=True)
    for g in sys.argv[1:] or GENERA:
        out = RAW / f"{g}.json"
        if out.exists():
            print(g, "cached"); continue
        try:
            recs = fetch_genus(g)
        except Exception as e:
            print(g, "FAILED", e); continue
        out.write_text(json.dumps(recs))
        print(g, len(recs), "records")
