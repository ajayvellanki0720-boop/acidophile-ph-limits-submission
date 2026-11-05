"""For each low-confidence Alicyclobacillaceae / Sulfolobaceae species, find the paper(s) behind its database value via open APIs
(Europe PMC): by DOI parsed from the ThermoBase source link where there is one, otherwise by species name (plus old names).
Publisher pages are not scraped (they serve bot challenges). Output: data/curated/primary_source_candidates.csv, for reading.
Usage: python3 scripts/18_primary_source_lookup.py"""
import html
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
tg = pd.read_csv(ROOT / "results/within_lineage_v2/primary_source_targets.csv")
tb = pd.read_csv(ROOT / "data/curated/thermobase_species.csv").set_index("species")
OLD = {"Saccharolobus solfataricus": ["Sulfolobus solfataricus"], "Saccharolobus shibatae": ["Sulfolobus shibatae"], "Saccharolobus caldissimus": ["Saccharolobus caldissimus"],
       "Sulfolobus tokodaii": ["Sulfolobus tokodaii"], "Alicyclobacillus cycloheptanicus": ["Bacillus cycloheptanicus", "Alicyclobacillus cycloheptanicus"],
       "Alicyclobacillus acidocaldarius": ["Bacillus acidocaldarius", "Alicyclobacillus acidocaldarius"], "Alicyclobacillus acidoterrestris": ["Bacillus acidoterrestris", "Alicyclobacillus acidoterrestris"]}
N = r"(\d+(?:\.\d+)?)"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + r"\s*(?:–|-|—|to|and)\s*" + N)


def epmc(q, n=15):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": n, "resultType": "core"})
    return json.load(urllib.request.urlopen(url, timeout=60))["resultList"]["result"]


def doi_from(url):
    m = re.search(r"(10\.\d{4,9}/[^\s#?&\"<>]+?)(?:#|\?|$|/-char|_pdf|\.pdf)", url + "#")
    return m.group(1) if m else None


rows = []
for r in tg.itertuples():
    sp = r.species
    srcs = str(tb.sources.get(sp, "")) if r.origin == "thermobase" else ""
    dois = [d for d in (doi_from(u.strip()) for u in srcs.split("|")) if d]
    pmcs = re.findall(r"PMC\d+", srcs)
    found = {}
    for d in dois:
        for x in epmc(f'DOI:"{d}"', 3):
            found[x.get("id")] = (x, f"ThermoBase source DOI {d}")
    for c in pmcs:
        for x in epmc(f'PMCID:{c}', 3):
            found[x.get("id")] = (x, f"ThermoBase source {c}")
    names = OLD.get(sp, [sp])
    for nm in names:
        q = f'(TITLE:"{nm}") AND (TITLE:"sp. nov" OR TITLE:"gen. nov" OR TITLE:"comb. nov" OR TITLE:"emended" OR TITLE:"description" OR TITLE:"new species" OR TITLE:"novel")'
        for x in epmc(q, 6):
            found.setdefault(x.get("id"), (x, f"title search '{nm}'"))
        time.sleep(0.2)
    for _id, (x, how) in found.items():
        ab = " ".join(html.unescape(re.sub(r"<[^>]+>", "", x.get("abstractText") or "")).split())
        sents = [s for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", ab) if "pH" in s]
        rng = next((m for s in sents for m in [RANGE.search(s)] if m), None)
        rows.append({"species": sp, "current_value": f"{r.ph_min}-{r.ph_max} (opt {r.ph_opt})", "via": how, "pmid": x.get("pmid"), "pmcid": x.get("pmcid"), "year": x.get("pubYear"),
                     "journal": x.get("journalTitle"), "title": re.sub(r"<[^>]+>", "", x.get("title") or "")[:120], "has_abstract": bool(ab),
                     "pH_sentence": (sents[0][:300] if sents else ""), "range_found": f"{rng.group(1)}-{rng.group(2)}" if rng else ""})
df = pd.DataFrame(rows).drop_duplicates(["species", "pmid", "title"])
df.to_csv(ROOT / "data/curated/primary_source_candidates.csv", index=False)
print(len(df), "candidate papers for", df.species.nunique(), "species;", int(df.has_abstract.sum()), "with abstracts;", int((df.range_found != "").sum()), "with a pH range in the abstract")
