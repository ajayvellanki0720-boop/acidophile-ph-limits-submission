"""Find a source for the 292 minima that have none, then classify them like the rest -> DATA_DIR/ph_min_recovered.csv
WHAT IS LEFT AFTER 107b. 2,710 of 3,002 minima cite a PMID and were checked against it. The remaining 292 cite a
database and nothing else: ThermoBase 217, AciDB 33, BacDive 15, and 27 parser rows whose PMID did not survive
merging. Those are the least documented rows in the study, and leaving them unchecked because the citation is
inconvenient would be choosing not to look.
TWO ROUTES, both open-API and neither touching a publisher's site:
  1. ThermoBase stores a source URL per species, and 88 of its 300 rows embed a DOI. Europe PMC resolves a DOI to
     its record, so those become checkable exactly like the PMID rows.
  2. For the rest there is no stored citation at all, so the original species description is SEARCHED for by name:
     Europe PMC, the binomial in quotes, restricted to the nomenclatural act that defines it ("sp. nov.").
     A hit counts only if the nomenclatural act in the title attaches to OUR binomial.
WHY THE TITLE RULE HAD TO BE TIGHTENED TWICE. Requiring only that the binomial appear somewhere in the title let two
distinct failure modes through, both found by auditing the 58 rows the first version recovered:
  1. Papers ABOUT the organism rather than describing it - "Thermostable lipases from ... Thermoanaerobacter
     thermohydrosulfuricus SOL1", whose pH range is an ENZYME activity range, and a plasmid-transformation paper on
     Acidiphilium multivorum. Fixed by requiring a nomenclatural marker (sp. nov., gen. nov., comb. nov.) in the title.
  2. RECLASSIFICATION papers, which that first fix does not catch, because they carry the marker and our binomial
     both - for a different reason. "Sulfuracidifex tepidarius gen. nov., sp. nov. and transfer of Sulfolobus
     metallicus ... to the genus Sulfuracidifex" names S. metallicus only as the taxon being moved; the pH range of
     0.4-5.5 in its abstract belongs to S. tepidarius. Likewise a Sulfurisphaera javensis description reclassifying
     Sulfolobus tokodaii. Reading either abstract as evidence about our species would attribute one organism's
     growth range to another, which is worse than recovering nothing.
     Fixed by requiring the binomial IMMEDIATELY PRECEDING the nomenclatural marker to be our species, so the paper
     has to be describing it rather than mentioning it.
WHAT A RECOVERED SOURCE DOES AND DOES NOT BUY. It lets the stored VALUE be checked against a paper, which is the
transcription question 107b answered for the other 2,710. It does not settle the BOUND, because the same genre
limit applies: these are abstracts, and abstracts state ranges rather than failed trials.
Usage: python3 scripts/107f_recover_missing_sources.py   (env DATA_DIR=results/broad2)"""
import json
import os
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
CACHE = ROOT / "data/raw/epmc_recover"
CACHE.mkdir(parents=True, exist_ok=True)
UA = "acidophile-ph-limits/1.0 (comparative genomics; open-API only)"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def epmc(query, n=4):
    key = re.sub(r"[^A-Za-z0-9]+", "_", query)[:120]
    fp = CACHE / f"{key}.json"
    if fp.exists():
        return json.loads(fp.read_text())
    url = f"{API}?{urllib.parse.urlencode({'query': query, 'resultType': 'core', 'format': 'json', 'pageSize': n})}"
    for a in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60) as r:
                j = json.load(r)
            break
        except Exception:
            if a == 2:
                j = {"resultList": {"result": []}}
            else:
                time.sleep(2 ** a)
    out = [{"pmid": x.get("pmid"), "doi": x.get("doi"), "title": x.get("title"),
            "abstract": x.get("abstractText"), "year": x.get("pubYear")}
           for x in j.get("resultList", {}).get("result", [])]
    fp.write_text(json.dumps(out, ensure_ascii=False))
    time.sleep(0.34)
    return out


ph = pd.read_csv(D / "phenotypes_all.csv")
ph["pmid"] = ph.source.astype(str).str.extract(r"PMID:(\d+)")
miss = ph[ph.pmid.isna()].copy()
print(f"1. {len(miss)} rows with no PMID")
print(miss.origin.value_counts().to_string())

tb = pd.read_csv(ROOT / "data/curated/thermobase_species.csv")
tb["doi"] = tb.sources.astype(str).str.extract(r"(10\.\d{4,9}/[^\s,;\"]+)")
doi_of = dict(zip(tb.species, tb.doi))
print(f"\n2. ThermoBase supplies a DOI for {tb.doi.notna().sum()} of its {len(tb)} species "
      f"({tb.doi.nunique()} distinct)")

NUM = r"([0-9]{1,2}(?:\.[0-9])?)"
RE_RANGE = re.compile(r"pH\s*(?:range\s*)?(?:of\s*|from\s*|between\s*)?" + NUM + r"\s*(?:-|to|and)\s*" + NUM, re.I)
RE_PHVAL = re.compile(r"pH[^.;]{0,40}?" + NUM, re.I)
RE_NOGROW = [
    re.compile(r"(?:no|not|without)\s+growth[^.;]{0,70}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"(?:did|was|were)\s+not\s+grow[^.;]{0,70}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
]


# The binomial each nomenclatural act in a title applies to: the capitalised-genus + epithet immediately before
# the marker. A title may contain several acts and several names, and only the name the act attaches to is the taxon
# the paper describes.
RE_ACT = re.compile(r"([A-Z][a-z]+)\s+([a-z]{3,})(?:\s+(?:subsp\.?\s+[a-z]+))?\s*,?\s*"
                    r"(?:gen\.?\s*nov\.?\s*,?\s*)?(?:sp\.?\s*nov|comb\.?\s*nov|subsp\.?\s*nov)")


def described(title):
    """Every binomial that a nomenclatural act in this title attaches to."""
    return [f"{m.group(1)} {m.group(2)}" for m in RE_ACT.finditer(title)]


def norm(t):
    t = unicodedata.normalize("NFKC", t or "")
    t = re.sub(r"<[^>]+>", " ", t)
    for a, b in (("–", "-"), ("—", "-"), ("−", "-"), (" ", " ")):
        t = t.replace(a, b)
    return re.sub(r"\s+", " ", t)


rows, n_doi, n_name, n_none = [], 0, 0, 0
for i, r in enumerate(miss.itertuples(), 1):
    sp = str(r.species)
    hit, route = None, ""
    d = doi_of.get(sp)
    if isinstance(d, str) and d:
        for h in epmc(f'DOI:"{d}"'):
            if h.get("abstract"):
                hit, route = h, "thermobase_doi"
                break
    if hit is None and len(sp.split()) >= 2:
        for h in epmc(f'"{sp}" AND "sp. nov."'):
            if not h.get("abstract"):
                continue
            if sp.lower() in {b.lower() for b in described(norm(h.get("title")))}:
                hit, route = h, "name_search"
                break
    if hit is None:
        n_none += 1
        rows.append({"species": sp, "ph_min": r.ph_min, "origin": r.origin, "route": "none",
                     "pmid": None, "verdict": "no_source_found", "evidence": ""})
        continue
    n_doi += route == "thermobase_doi"
    n_name += route == "name_search"
    t = norm((hit.get("title") or "") + ". " + (hit.get("abstract") or ""))
    nogrow = [(float(m.group(1)), m.start()) for rx in RE_NOGROW for m in rx.finditer(t)]
    ranges = [(float(m.group(1)), float(m.group(2)), m.start()) for m in RE_RANGE.finditer(t)
              if float(m.group(1)) < float(m.group(2)) <= 14]
    lo = [g for g in ranges if abs(g[0] - r.ph_min) <= 0.25]
    val = [m for m in RE_PHVAL.finditer(t) if abs(float(m.group(1)) - r.ph_min) <= 0.05]
    below = [(v, p) for v, p in nogrow if v <= r.ph_min + 1e-9]

    def s(pos):
        a = max(t.rfind(". ", 0, pos) + 1, 0); b = t.find(". ", pos)
        return t[a:(b + 1 if b != -1 else len(t))].strip()[:350]

    if below:
        v, p = max(below); verdict, ev = "verified_measured", s(p)
    elif lo:
        verdict, ev = "range_only", s(lo[0][2])
    elif val:
        verdict, ev = "value_stated", s(val[0].start())
    elif ranges:
        verdict, ev = "value_not_in_source", s(ranges[0][2])
    else:
        verdict, ev = "no_ph_statement", ""
    rows.append({"species": sp, "ph_min": r.ph_min, "origin": r.origin, "route": route,
                 "pmid": hit.get("pmid"), "verdict": verdict, "evidence": ev})
    if i % 50 == 0:
        print(f"   {i}/{len(miss)} processed", flush=True)

rec = pd.DataFrame(rows)
print(f"\n2b. title rule: a hit is accepted only when a nomenclatural act in its title attaches to the species "
      f"itself, which rejects papers that merely study the organism and reclassification papers that move it")
print(f"\n3. source recovered for {n_doi + n_name} of {len(miss)} "
      f"({(n_doi + n_name) / len(miss):.1%}): {n_doi} by ThermoBase DOI, {n_name} by name search; "
      f"{n_none} still have none")
print("\n4. verdicts on the recovered rows")
for k, n in rec.verdict.value_counts().items():
    print(f"   {k:22s} {n:4d}  {n / len(rec):6.1%}")
bad = rec[rec.verdict == "value_not_in_source"]
print(f"\n5. {len(bad)} recovered rows state a pH range that does not contain the stored minimum")
for r in bad.head(12).itertuples():
    print(f"   {str(r.species)[:38]:38s} stored={r.ph_min}")
    print(f"      \"{str(r.evidence)[:150]}\"")
rec.to_csv(D / "ph_min_recovered.csv", index=False)
print(f"\nwrote {D / 'ph_min_recovered.csv'}")
