"""Fetch species-description abstracts from Europe PMC into data/raw/europepmc/<key>.json (raw, never edited).
Usage: python3 scripts/00c_fetch_europepmc.py"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "europepmc"
GENERA = """Acidithiobacillus Leptospirillum Acidiphilium Acidocella Acidisphaera Acidisoma Acidicaldus Acidobacterium
Granulicella Terriglobus Edaphobacter Acidipila Bryocella Sulfobacillus Alicyclobacillus Acidimicrobium Ferrimicrobium
Ferrithrix Acidithrix Aciditerrimonas Acidothermus Thiomonas Acidiferrobacter Acidihalobacter Acidibacillus Acidomonas
Sulfolobus Saccharolobus Acidianus Metallosphaera Sulfurisphaera Stygiolobus Sulfuracidifex Picrophilus Thermoplasma
Ferroplasma Cuniculiplasma Acidiplasma Thermogymnomonas Aciduliprofundum Acidilobus Caldisphaera Vulcanisaeta
Acidiella Acidocaldus Alicyclobacter Hydrogenobaculum Methylacidiphilum Methylacidimicrobium Acidimicrobiaceae
Desulfosporosinus Acidiphilum Tepidimonas Mycobacterium_acidophilic Bacillus_acidocaldarius Oleispira Gluconacetobacter
Acetobacter Gluconobacter Lactobacillus_acid Bifidobacterium_acid Helicobacter Streptococcus_mutans Escherichia_acid
Thermoanaerobacter_acid Acidovorax Rhodopila Acidiphilium_like Sulfuriferula Acidithrix_like Ferroacidibacillus Tumebacter""".split()
GENERA = [g for g in GENERA if "_" not in g and not g.endswith("like")]
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch(query, cap=1000):
    out, cursor = [], "*"
    while len(out) < cap:
        url = BASE + "?" + urllib.parse.urlencode({"query": query, "format": "json", "pageSize": 100,
                                                    "resultType": "core", "cursorMark": cursor})
        d = json.load(urllib.request.urlopen(url, timeout=60))
        res = d["resultList"]["result"]
        out += [{k: r.get(k) for k in ("id", "pmid", "doi", "title", "abstractText", "journalTitle", "pubYear", "authorString")} for r in res]
        nxt = d.get("nextCursorMark")
        if not res or not nxt or nxt == cursor:
            break
        cursor = nxt
        time.sleep(0.3)
    return out


if __name__ == "__main__":
    RAW.mkdir(parents=True, exist_ok=True)
    qs = {g: f'ABSTRACT:"{g}" AND ABSTRACT:"pH" AND (ABSTRACT:"sp. nov." OR ABSTRACT:"gen. nov." OR ABSTRACT:"emended description" OR ABSTRACT:"comb. nov.")' for g in GENERA}
    qs["_general"] = '(ABSTRACT:"acidophilic" OR ABSTRACT:"acidophile" OR ABSTRACT:"extremely acidophilic") AND ABSTRACT:"pH" AND (ABSTRACT:"sp. nov." OR ABSTRACT:"gen. nov.")'
    for k, q in qs.items():
        f = RAW / f"{k}.json"
        if f.exists():
            continue
        try:
            recs = fetch(q)
        except Exception as e:
            print(k, "FAILED", e); continue
        f.write_text(json.dumps(recs)); print(k, len(recs), flush=True)
