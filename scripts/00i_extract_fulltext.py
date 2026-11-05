"""Extract growth-pH statements and explicit-bound statements from cached full texts -> data/curated/fulltext_review.csv.
One row per paper with: all growth-pH range sentences (short snippets), and bound-evidence sentences such as
'no growth at pH 1.5' / 'did not grow below pH 2'. Human reading assigns species, values and bound type afterwards."""
import glob
import html
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N = r"(\d+(?:\.\d+)?)"
SEP = r"\s*(?:–|-|—|−|to|and)\s*"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + SEP + N)
BOUND = re.compile(r"(no|not|never|failed to|did not|didn't|unable to|without|absence of|lack of)[^.;]{0,50}(grow|growth)[^.;]{0,60}pH|"
                   r"(grow|growth)[^.;]{0,40}(was|were|is)?\s*(not|never)\s*(observed|detected|occur)[^.;]{0,50}pH|"
                   r"pH[^.;]{0,40}(no growth|not grow|did not|failed)|(below|above|lower than|higher than)\s*pH\s*" + N, re.I)
GROW = re.compile(r"grow|growth|grew", re.I)


def text_of(path):
    root = ET.parse(path).getroot()
    title = " ".join("".join(root.find(".//article-title").itertext()).split()) if root.find(".//article-title") is not None else ""
    paras = ["".join(p.itertext()) for p in root.iter("p")]
    return title, [" ".join(html.unescape(p).split()) for p in paras]


rows = []
listing = {i["pmcid"]: i for i in json.load(open(ROOT / "data/curated/fulltext_listing.json"))}
for f in sorted(glob.glob(str(ROOT / "data/raw/europepmc_ft/*.xml"))):
    pmcid = Path(f).stem
    try:
        title, paras = text_of(f)
    except ET.ParseError:
        continue
    ranges, bounds = [], []
    for p in paras:
        for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", p):
            if "pH" not in s or not GROW.search(s):
                continue
            m = RANGE.search(s)
            if m and float(m.group(1)) < float(m.group(2)):
                ranges.append(s[max(0, m.start() - 70):m.start() + 170])
            if BOUND.search(s):
                bounds.append(s[:230])
    if not ranges and not bounds:
        continue
    it = listing.get(pmcid, {})
    rows.append({"pmcid": pmcid, "pmid": it.get("pmid"), "year": it.get("pubYear"), "title": title[:140],
                 "n_range": len(ranges), "n_bound": len(bounds), "ranges": " ## ".join(ranges[:4]), "bounds": " ## ".join(bounds[:4])})
out = pd.DataFrame(rows)
out.to_csv(ROOT / "data/curated/fulltext_review.csv", index=False)
print(len(out), "papers with growth-pH statements |", (out.n_bound > 0).sum(), "with explicit-bound wording")
