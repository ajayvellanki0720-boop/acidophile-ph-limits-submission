"""Targeted search for well-studied acidophiles missing from species-description abstracts.
For each species, pull Europe PMC abstracts that name it and contain a pH range with growth wording;
output data/curated/target_species_review.csv for manual reading (nothing is auto-accepted)."""
import html
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPECIES = """Picrophilus torridus|Picrophilus oshimae|Thermoplasma acidophilum|Thermoplasma volcanium|Ferroplasma acidarmanus|
Acidianus brierleyi|Acidianus ambivalens|Acidianus infernus|Metallosphaera sedula|Metallosphaera yellowstonensis|Metallosphaera cuprina|
Sulfolobus acidocaldarius|Sulfolobus solfataricus|Saccharolobus solfataricus|Sulfolobus tokodaii|Sulfolobus islandicus|Acidilobus aceticus|
Leptospirillum ferrooxidans|Leptospirillum ferriphilum|Leptospirillum ferrodiazotrophum|Acidithiobacillus ferrooxidans|Acidithiobacillus thiooxidans|
Acidithiobacillus caldus|Acidithiobacillus ferrivorans|Acidiphilium cryptum|Acidiphilium multivorum|Acidiphilium acidophilum|Acidimicrobium ferrooxidans|
Sulfobacillus acidophilus|Sulfobacillus thermosulfidooxidans|Alicyclobacillus acidocaldarius|Alicyclobacillus acidoterrestris|Acidocella facilis|
Acidobacterium capsulatum|Methylacidiphilum infernorum|Methylacidiphilum fumariolicum|Hydrogenobaculum acidophilum|Acidisphaera rubrifaciens|
Acidiferrobacter thiooxydans|Ferrimicrobium acidiphilum|Ferroplasma acidiphilum|Halorubrum|Acidihalobacter prosperus|Thermoplasma|Cuniculiplasma divulgatum|
Sulfurisphaera ohwakuensis|Stygiolobus azoricus|Metallosphaera hakonensis|Acidianus copahuensis|Sulfolobus metallicus|Vulcanisaeta distributa""".replace("\n", "").split("|")
N = r"(\d+(?:\.\d+)?)"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + r"\s*(?:–|-|—|to|and)\s*" + N)
rows = []
for sp in SPECIES:
    if not sp or sp == "Halorubrum":
        continue
    q = f'"{sp}" AND (pH) AND (growth OR grows OR grew)'
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 60, "resultType": "core"})
    try:
        res = json.load(urllib.request.urlopen(url, timeout=60))["resultList"]["result"]
    except Exception as e:
        print(sp, "FAILED", e); continue
    for r in res:
        ab = " ".join(html.unescape(re.sub(r"<[^>]+>", "", r.get("abstractText") or "")).split())
        if sp.split()[-1] not in ab and sp.split()[0] not in ab:
            continue
        for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", ab):
            m = RANGE.search(s)
            if m and float(m.group(1)) < float(m.group(2)) and re.search(r"grow|growth|grew", s, re.I):
                rows.append({"species": sp, "pmid": r.get("pmid"), "year": r.get("pubYear"), "title": (r.get("title") or "")[:90], "lo": m.group(1), "hi": m.group(2), "sentence": s[:260]})
                break
pd.DataFrame(rows).to_csv(ROOT / "data/curated/target_species_review.csv", index=False)
print(len(rows), "candidate statements for", len({r["species"] for r in rows}), "species")
