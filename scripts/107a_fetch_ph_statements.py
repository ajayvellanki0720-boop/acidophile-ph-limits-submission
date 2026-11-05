"""Fetch the source abstract for every reported minimum that carries a PMID -> data/raw/epmc_verify/*.json
WHY THIS EXISTS. docs/MANUSCRIPT.md names the outcome variable as the study's principal limitation, and the specific
defect is narrow and known: `ph_min_bound` reads "measured" for 2,732 of 3,002 genomes because that string is the
DEFAULT the parser writes (scripts/00g line 72), not something read from the paper. 97% of those labels were never
checked against a source. A label nobody read is not evidence, and the honest fix is to go back to the papers.
WHAT IT FETCHES. Europe PMC's REST search endpoint, which is open and needs no key, queried in batches of 40 PMIDs
with resultType=core so each hit carries abstractText. 2,708 distinct PMIDs back 2,710 phenotype rows. Results are
cached per batch so a re-run costs nothing and the classification in 107b can be re-run freely as its rules change.
No publisher page is touched and no bot challenge is answered: abstracts of species descriptions carry the growth
range, which is the sentence this needs, and where they do not the row is reported unresolved rather than guessed.
Usage: python3 scripts/107a_fetch_ph_statements.py   (env DATA_DIR=results/broad2)"""
import json
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
CACHE = ROOT / "data/raw/epmc_verify"
CACHE.mkdir(parents=True, exist_ok=True)
UA = "acidophile-ph-limits/1.0 (comparative genomics; open-API use only)"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

ph = pd.read_csv(D / "phenotypes_all.csv")
ph["pmid"] = ph.source.astype(str).str.extract(r"PMID:(\d+)")
pmids = sorted(ph.pmid.dropna().unique())
print(f"1. {len(ph)} phenotype rows; {ph.pmid.notna().sum()} carry a PMID; {len(pmids)} distinct PMIDs")

BATCH = 40
batches = [pmids[i:i + BATCH] for i in range(0, len(pmids), BATCH)]
print(f"2. fetching in {len(batches)} batches of up to {BATCH}")
n_new = n_cached = 0
for bi, b in enumerate(batches):
    fp = CACHE / f"batch_{bi:04d}.json"
    if fp.exists():
        n_cached += 1
        continue
    q = " OR ".join(f"EXT_ID:{p}" for p in b)
    url = f"{API}?{urllib.parse.urlencode({'query': q, 'resultType': 'core', 'format': 'json', 'pageSize': BATCH + 10})}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=90) as r:
                payload = json.load(r)
            break
        except Exception as e:
            if attempt == 3:
                raise
            print(f"   batch {bi}: {type(e).__name__}, retrying in {2 ** attempt}s")
            time.sleep(2 ** attempt)
    # Only the fields the classifier reads are kept. Storing whole records would make the cache enormous and would
    # also store full texts whose licences differ; abstract + title + identifiers is what a growth range needs.
    keep = [{"pmid": r.get("pmid"), "pmcid": r.get("pmcid"), "doi": r.get("doi"), "title": r.get("title"),
             "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"), "year": r.get("pubYear"),
             "isOpenAccess": r.get("isOpenAccess"), "abstract": r.get("abstractText")}
            for r in payload.get("resultList", {}).get("result", [])]
    fp.write_text(json.dumps({"requested": b, "result": keep}, ensure_ascii=False))
    n_new += 1
    if n_new % 10 == 0:
        print(f"   {n_new} fetched, {bi + 1}/{len(batches)} batches done")
    time.sleep(0.34)          # Europe PMC asks for courteous use; this is ~3 req/s at the batch level
print(f"3. {n_new} batches fetched, {n_cached} already cached")

got, abstracts = {}, 0
for fp in sorted(CACHE.glob("batch_*.json")):
    for r in json.loads(fp.read_text())["result"]:
        if r.get("pmid"):
            got[r["pmid"]] = r
            abstracts += bool(r.get("abstract"))
missing = [p for p in pmids if p not in got]
print(f"4. {len(got)} of {len(pmids)} PMIDs returned a record; {abstracts} of those carry an abstract")
print(f"   {len(missing)} returned nothing from Europe PMC")
# A row whose source yields no text cannot be verified either way, and that is a reportable outcome rather than a
# failure to route around. It is recorded so 107b counts it as unresolved instead of silently dropping it.
pd.DataFrame({"pmid": missing}).to_csv(D / "ph_min_sources_missing.csv", index=False)
print(f"wrote {D / 'ph_min_sources_missing.csv'}")
