"""Verify the bound from OPEN-ACCESS FULL TEXT, and measure the ceiling on doing so -> DATA_DIR/ph_min_fulltext.csv
WHY ABSTRACTS WERE NOT ENOUGH, which 107b established rather than assumed. Of 2,708 source abstracts, 94.9% state a
pH range and 84.9% an optimum, but only 1.1% contain the phrase "no growth" at all and NONE contains "lowest pH
tested". Inspection of the no-growth sentences that do occur shows they are almost all about temperature, salt or
carbon source. A species-description abstract states the range and leaves the bracketing to the paper's tables, so
an abstract is structurally incapable of settling whether a minimum is measured or censored. That is a property of
the genre, not a failure of the regexes.
THE CEILING, stated before the result because it is the main finding. Europe PMC reports open-access full text for
68 of the 2,708 sources - 2.5%. The remainder are in subscription journals, overwhelmingly IJSEM and its
predecessors, where the species descriptions of most prokaryotes live. The project's standing rule is to use open
APIs only and not to defeat publisher access controls, so 2.5% is the honest ceiling for full-text verification by
this route, and nothing here should be read as a sample of the other 97.5%: open-access status correlates with
journal, era and funder, so the verified subset is not representative.
WHAT IT DOES. Fetches the OA full text for those sources from Europe PMC's textMinedTerms/fullTextXML endpoint,
searches the whole article rather than the abstract, and applies the same precedence as 107b with one addition: a
full text can state the tested grid ("pH 2.0, 3.0, 4.0 ... were tested"), which settles censoring directly, so a
grid whose minimum equals the reported minimum is explicit censoring evidence.
Usage: python3 scripts/107c_fulltext_bound.py   (env DATA_DIR=results/broad2)"""
import json
import os
import re
import time
import unicodedata
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
FT = ROOT / "data/raw/epmc_fulltext"
FT.mkdir(parents=True, exist_ok=True)
UA = "acidophile-ph-limits/1.0 (comparative genomics; open-API use only)"

oa = []
for fp in sorted((ROOT / "data/raw/epmc_verify").glob("batch_*.json")):
    for r in json.loads(fp.read_text())["result"]:
        if r.get("pmid") and r.get("pmcid") and r.get("isOpenAccess") == "Y":
            oa.append((r["pmid"], r["pmcid"]))
print(f"1. {len(oa)} sources with open-access full text (of 2,708 = {len(oa) / 2708:.1%})")

n_new = n_cached = n_fail = 0
for pmid, pmcid in oa:
    fp = FT / f"{pmcid}.xml"
    if fp.exists() and fp.stat().st_size > 500:
        n_cached += 1
        continue
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=90) as r:
            fp.write_bytes(r.read())
        n_new += 1
    except Exception as e:
        # A refusal here is a licence or availability fact about that article, not something to route around.
        n_fail += 1
        print(f"   {pmcid}: {type(e).__name__} - left unfetched")
    time.sleep(0.34)
print(f"2. {n_new} fetched, {n_cached} cached, {n_fail} unavailable")

# Full texts gathered earlier in the project for a different purpose are reused where the article matches.
legacy = {p.stem: p for p in (ROOT / "data/raw/europepmc_ft").glob("PMC*.xml")}
print(f"3. {len(legacy)} full texts already in data/raw/europepmc_ft; "
      f"{len({p for _, p in oa} & set(legacy))} of them are among the sources needed here")


def norm(t):
    t = unicodedata.normalize("NFKC", t or "")
    t = re.sub(r"<[^>]+>", " ", t)
    for a, b in (("–", "-"), ("—", "-"), ("−", "-"), ("‐", "-"),
                 ("≥", ">="), ("≤", "<="), (" ", " ")):
        t = t.replace(a, b)
    return re.sub(r"\s+", " ", t)


NUM = r"([0-9]{1,2}(?:\.[0-9])?)"
RE_NOGROW = [
    re.compile(r"(?:no|not|without)\s+growth[^.;]{0,90}?(?:at|below|under|beneath)\s*(?:a\s*)?pH\s*(?:of\s*|value[s]?\s*(?:of\s*)?)?" + NUM, re.I),
    re.compile(r"(?:did|does|could|was|were)\s+not\s+grow[^.;]{0,90}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"growth\s+(?:was|were)?\s*(?:not|never)\s+(?:observed|detected|found|evident)[^.;]{0,90}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"pH\s*(?:of\s*)?" + NUM + r"[^.;]{0,60}?(?:no|not)\s+(?:growth|grow)", re.I),
]
RE_CENSOR = [
    re.compile(r"lowest\s+(?:pH\s+)?(?:value\s+)?(?:level\s+)?(?:tested|examined|assayed|investigated|studied)", re.I),
    re.compile(r"(?:down\s+to|as\s+low\s+as)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"pH\s*(?:>=|at\s+least)\s*" + NUM, re.I),
]
# A stated tested grid: three or more pH values enumerated as what was assayed. If its minimum equals the reported
# minimum, nothing below the limit was tried and the value is censored by construction.
RE_GRID = re.compile(r"pH\s*(?:values?\s*)?(?:of\s*)?((?:[0-9]{1,2}(?:\.[0-9])?\s*,\s*){2,}[0-9]{1,2}(?:\.[0-9])?)[^.;]{0,60}?(?:test|examin|assay|investigat|incubat|grown)", re.I)

ph = pd.read_csv(D / "phenotypes_all.csv")
ph["pmid"] = ph.source.astype(str).str.extract(r"PMID:(\d+)")
pm2pmc = dict(oa)

rows = []
for r in ph[ph.pmid.isin(pm2pmc)].itertuples():
    pmcid = pm2pmc[r.pmid]
    fp = FT / f"{pmcid}.xml"
    if not (fp.exists() and fp.stat().st_size > 500):
        fp = legacy.get(pmcid)
    if fp is None or not fp.exists():
        rows.append({"species": r.species, "pmid": r.pmid, "pmcid": pmcid, "ph_min": r.ph_min,
                     "verdict": "fulltext_unavailable", "evidence": ""})
        continue
    t = norm(fp.read_text(errors="replace"))
    nogrow = [(float(m.group(1)), m.start()) for rx in RE_NOGROW for m in rx.finditer(t)]
    censor = [(m.group(0), m.start()) for rx in RE_CENSOR for m in rx.finditer(t)]
    grids = []
    for m in RE_GRID.finditer(t):
        vals = [float(x) for x in re.findall(r"[0-9]{1,2}(?:\.[0-9])?", m.group(1))]
        if vals and max(vals) <= 14:
            grids.append((min(vals), m.start()))
    below = [(v, p) for v, p in nogrow if v <= r.ph_min + 1e-9]

    def s(pos):
        a = max(t.rfind(". ", 0, pos) + 1, 0); b = t.find(". ", pos)
        return t[a:(b + 1 if b != -1 else len(t))].strip()[:400]

    if below:
        v, p = max(below)
        rows.append({"species": r.species, "pmid": r.pmid, "pmcid": pmcid, "ph_min": r.ph_min,
                     "verdict": "verified_measured", "evidence": s(p), "nogrow_ph": v})
    elif [g for g in grids if abs(g[0] - r.ph_min) <= 0.25]:
        g = [g for g in grids if abs(g[0] - r.ph_min) <= 0.25][0]
        rows.append({"species": r.species, "pmid": r.pmid, "pmcid": pmcid, "ph_min": r.ph_min,
                     "verdict": "verified_censored", "evidence": s(g[1])})
    elif censor:
        rows.append({"species": r.species, "pmid": r.pmid, "pmcid": pmcid, "ph_min": r.ph_min,
                     "verdict": "verified_censored", "evidence": s(censor[0][1])})
    else:
        rows.append({"species": r.species, "pmid": r.pmid, "pmcid": pmcid, "ph_min": r.ph_min,
                     "verdict": "fulltext_silent", "evidence": ""})

f = pd.DataFrame(rows)
print(f"\n4. verdicts on the {len(f)} phenotype rows whose source has open full text")
for k, n in f.verdict.value_counts().items():
    print(f"   {k:24s} {n:4d}  {n / len(f):6.1%}")
f.to_csv(D / "ph_min_fulltext.csv", index=False)
print(f"\nwrote {D / 'ph_min_fulltext.csv'}")
