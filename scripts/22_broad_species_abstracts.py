"""Broad Europe PMC sweep: every species/genus description abstract (2005 onward) that mentions a pH range, regardless of acid-related
wording. Keeps compact records in data/raw/europepmc_broad/<year>.json (git-ignored). Usage: python3 scripts/22_broad_species_abstracts.py [first_year] [last_year]"""
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/europepmc_broad"; OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
N = r"(\d+(?:\.\d+)?)"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + r"\s*(?:–|-|—|−|to|and)\s*" + N)
DESC = re.compile(r"sp\. nov|gen\. nov|comb\. nov|nov\.", re.I)


def get(url):
    for i in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=90))
        except Exception:
            time.sleep(3 * (i + 1))
    raise RuntimeError(url)


y0, y1 = int(sys.argv[1]) if len(sys.argv) > 1 else 2005, int(sys.argv[2]) if len(sys.argv) > 2 else 2026
for y in range(y0, y1 + 1):
    f = OUT / f"{y}.json"
    if f.exists():
        continue
    q = f'("sp. nov." OR "gen. nov." OR "comb. nov.") AND pH AND PUB_YEAR:{y} AND (growth OR grew OR grows)'
    cursor, kept, seen = "*", [], 0
    while True:
        d = get(BASE + "?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 100, "resultType": "core", "cursorMark": cursor}))
        res = d["resultList"]["result"]
        for r in res:
            seen += 1
            title = html.unescape(re.sub(r"<[^>]+>", "", r.get("title") or ""))
            ab = html.unescape(re.sub(r"<[^>]+>", "", r.get("abstractText") or ""))
            if DESC.search(title) and RANGE.search(ab):
                kept.append({"id": r.get("id"), "pmid": r.get("pmid"), "doi": r.get("doi"), "pmcid": r.get("pmcid"), "year": r.get("pubYear"),
                             "journal": r.get("journalTitle"), "title": " ".join(title.split()), "abstract": " ".join(ab.split())})
        nxt = d.get("nextCursorMark")
        if not res or not nxt or nxt == cursor:
            break
        cursor = nxt
        time.sleep(0.25)
    f.write_text(json.dumps(kept))
    print(y, f"{seen} abstracts scanned, {len(kept)} species descriptions with a pH range", flush=True)
