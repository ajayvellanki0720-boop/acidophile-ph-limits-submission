"""Can a literature index be probed for demonstrated growth limits? No, and this records why -> DATA_DIR/demonstrated_limits.csv
The project's largest remaining weakness is that only 5 species have a demonstrated limit in the strict sense: a reported minimum is usually the lowest
pH at which growth was seen, not a pH below which growth was shown to fail. Upgrading that needs a statement of failure to grow, and most species
descriptions are paywalled, so the idea tested here was to probe Europe PMC's INDEX rather than fetch any text - ask whether a given paper contains
"no growth at pH 4.0" without ever reading it, which costs nothing and touches no paywall. Batched 50 PMIDs at a time it screened all 2,708 papers
behind the phenotype table in a few hundred queries and appeared to work: 205 papers (7.6%) came back positive, a value was recovered for 98 of them,
and the subset gave an estimate of +0.52 per SD against the headline's +0.185.
ALL OF THAT WAS AN ARTEFACT AND NONE OF IT IS REPORTED, because the index strips negations as stopwords. The two checks below are what established it,
and they are kept as the script's whole content so that nobody - including a later me - spends another afternoon on this approach.
  1. THE MECHANISM. `"no growth at pH 4.0"`, `"growth at pH 4.0"` and `"not growth at pH 4.0"` return the identical hit count, while `"never growth at
     pH 4.0"` and a nonsense word in the same slot return zero. So "no" and "not" are dropped from the phrase and every negated query silently matches
     the affirmative statement. A query for papers reporting a failure to grow returns papers reporting growth.
  2. GROUND TRUTH HELD LOCALLY. The project already stores the statement each parser row was built from, so the probe can be checked without any
     network at all. Of the papers the screen called positive and whose statement is stored, almost none contain a failure-to-grow clause of any kind;
     the "recovered" pH is simply the recorded minimum read back out of a sentence like "capable of growth at pH 4.0-7.2".
That also explains the inflated estimate, which should have been the warning sign: the positive set was papers whose text states a growth bound at one
of the probed values, and since the probed grid stopped at pH 6 it was enriched for acidophiles - the same enrichment confound that the evidence-tier
test (scripts/58) had already shown produces estimates two to three times the headline. An estimate that large from a subset selected on the outcome's
own range should have prompted the validation before the write-up rather than after.
No earlier script in this project used a negated index query - checked across every script that calls the Europe PMC API - so the hazard is confined to
this one and was caught before anything was published.
WHAT WOULD ACTUALLY WORK. Only reading the sentence settles it, which means full text. Europe PMC serves full text for a small minority of species
descriptions and the project already harvests those (`scripts/00h`, `scripts/00i`), where the local parse does find genuine bounds because it reads the
negation. Extending the count therefore needs either more open full text than exists or hand-reading of paywalled descriptions, which is where the
synthesis already places this gap. It stays unquantified, and that is the honest position.
Usage: python3 scripts/65_demonstrated_limits.py   (env DATA_DIR=results/broad2)"""
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
CACHE = ROOT / "data/raw/epmc_phrase"
CACHE.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
rows = []


def hits(query):
    """Hit count for one Europe PMC search, cached on disk by query hash."""
    key = CACHE / ("n_" + hashlib.sha1(query.encode()).hexdigest() + ".json")
    if key.exists():
        return json.loads(key.read_text())
    url = f"{BASE}/search?query={urllib.parse.quote(query)}&format=json&pageSize=1&resultType=idlist"
    for attempt in range(4):
        try:
            n = json.load(urllib.request.urlopen(url, timeout=120))["hitCount"]
            break
        except Exception as e:
            if attempt == 3:
                print("   (query failed:", e, ")")
                return None
            time.sleep(2 ** attempt)
    key.write_text(json.dumps(n))
    time.sleep(0.25)
    return n


print("1. the index strips negations, so a negated phrase query matches the affirmative statement")
PAIRS = ['"no growth at pH 4.0"', '"growth at pH 4.0"', '"not growth at pH 4.0"',
         '"never growth at pH 4.0"', '"xylophone growth at pH 4.0"',
         '"no growth at pH"', '"growth at pH"']
counts = {}
for q in PAIRS:
    counts[q] = hits(q)
    print(f"   {q:36s} {counts[q]}")
    rows.append({"part": "stopword test", "query": q, "hit_count": counts[q]})
same = counts.get('"no growth at pH 4.0"') == counts.get('"growth at pH 4.0"') == counts.get('"not growth at pH 4.0"')
print(f"   negated and affirmative counts identical: {same}; a non-stopword word in the same slot returns "
      f"{counts.get(chr(34) + 'never growth at pH 4.0' + chr(34))}")
# the same thing on one paper whose stored statement reports growth at pH 4.0 and no failure at all
for q in ['EXT_ID:24867171 AND "no growth at pH 4.0"', 'EXT_ID:24867171 AND "growth at pH 4.0"',
          'EXT_ID:24867171 AND "no growth at pH 9.9"']:
    n = hits(q)
    print(f"   {q:46s} {n}")
    rows.append({"part": "stopword test", "query": q, "hit_count": n})
print("   That paper's stored statement reads \"capable of growth at pH 4.0-7.2 (optimum pH 5.5-6.0)\" and contains no negation,")
print("   yet it answers a query for a failure to grow at pH 4.0.")

print("\n2. the same conclusion from ground truth the project already holds, with no network at all")
cand = pd.read_csv(ROOT / "data/curated/broad_candidates.csv")
cand["pmid"] = cand.pmid.astype(str)
stored = cand.dropna(subset=["statement"]).groupby("pmid").statement.apply(lambda s: " ## ".join(s.astype(str)))
# the 205 papers the index screen called positive, recovered from the cached screen results
screened = set()
for f in CACHE.glob("*.json"):
    if f.name.startswith("n_"):
        continue
    try:
        screened |= set(json.loads(f.read_text()))
    except Exception:
        pass
NOGROW = re.compile(r"(no growth|did not grow|not observed|failed to grow|unable to grow)[^.]{0,60}?pH\s*(\d+(?:\.\d+)?)", re.I)
check = stored[stored.index.isin(screened)]
with_clause = int(check.map(lambda s: bool(NOGROW.search(s))).sum())
print(f"   {len(check)} of the screened-positive papers have their statement stored in this repository")
print(f"   of those, {with_clause} contain a failure-to-grow clause of any kind ({with_clause / max(len(check), 1):.1%})")
print("   so the screen's positives are papers that state a growth RANGE at a probed value, not a failure to grow")
rows.append({"part": "ground truth", "n_checked": len(check), "n_with_nogrowth_clause": with_clause,
             "rate": with_clause / max(len(check), 1)})

print("\n3. for comparison, the same regex on full text the project did fetch and parse locally")
fr = ROOT / "data/curated/fulltext_review.csv"
if fr.exists():
    ft = pd.read_csv(fr)
    n_bound = int((ft.n_bound.fillna(0) > 0).sum())
    print(f"   {len(ft)} open-access papers with full text fetched by scripts/00h; {n_bound} carry a bound statement the local parse found")
    print("   Reading the text works because the parse sees the negation. There are 99 such papers, against 2,708 in the dataset.")
    rows.append({"part": "fulltext comparison", "n_papers": len(ft), "n_with_bound": n_bound})

pd.DataFrame(rows).to_csv(D / "demonstrated_limits.csv", index=False)
print("\nwrote", D / "demonstrated_limits.csv")
print("\nConclusion. The index cannot be asked this question. The count of demonstrated limits stays at the 5 that were read by hand,")
print("the weakness stays unquantified, and the synthesis's entry on it is unchanged except to record that this route is closed.")
print("The transferable lesson is narrower than 'validate your data source': a full-text search engine is not a proposition checker,")
print("and a query whose meaning depends on a function word cannot be trusted to a search index at all. The warning sign was present")
print("before the validation - a subset selected on the outcome's own range returned an estimate 2.8x the headline, which this project")
print("had already learned to read as enrichment rather than signal.")
