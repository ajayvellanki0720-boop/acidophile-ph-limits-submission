"""Fetch Europe PMC full-text XML for acid-related species descriptions (raw -> data/raw/europepmc_ft/<PMCID>.xml).
Only what Europe PMC serves is stored; rows without retrievable text are logged. Raw files are git-ignored."""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/europepmc_ft"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
TERMS = ['(acidophilic OR acidophile OR acidotolerant OR acidiphilic OR thermoacidophilic OR "acid mine drainage" OR "acidic hot spring" OR "acidic soil" OR "acidic peat")',
         '("sp. nov." OR "gen. nov." OR "comb. nov.")', "HAS_FT:y"]


def listing():
    out, cursor = [], "*"
    q = " AND ".join(TERMS)
    while True:
        url = BASE + "/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 100, "resultType": "lite", "cursorMark": cursor})
        d = json.load(urllib.request.urlopen(url, timeout=60))
        res = d["resultList"]["result"]
        out += [{k: r.get(k) for k in ("id", "pmid", "pmcid", "doi", "title", "pubYear", "isOpenAccess")} for r in res]
        nxt = d.get("nextCursorMark")
        if not res or not nxt or nxt == cursor:
            return out
        cursor = nxt


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    items = [i for i in listing() if i.get("pmcid")]
    (ROOT / "data/curated/fulltext_listing.json").write_text(json.dumps(items))
    ok = bad = 0
    for it in items:
        f = OUT / f"{it['pmcid']}.xml"
        if f.exists():
            ok += 1; continue
        try:
            f.write_bytes(urllib.request.urlopen(f"{BASE}/{it['pmcid']}/fullTextXML", timeout=60).read()); ok += 1
        except Exception:
            bad += 1
        time.sleep(0.25)
    print(f"{len(items)} listed with PMCID | {ok} full texts stored | {bad} unavailable")
