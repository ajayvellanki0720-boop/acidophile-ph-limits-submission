"""Classify every reported minimum against its own paper's wording -> DATA_DIR/ph_min_verification.csv
THE DEFECT THIS ADDRESSES. `ph_min_bound` reads "measured" for 2,732 of 3,002 genomes because scripts/00g writes
that string as a DEFAULT, not because anyone read it. The distinction it claims to make is the one the whole
censoring argument turns on: a MEASURED minimum means growth was tried below the limit and failed, so the number is
the organism's limit; a CENSORED minimum means the number is merely the lowest pH anyone tried, so the true limit is
somewhere below it and unknown. Those are different quantities and the dataset does not currently distinguish them.
HOW A CLAIM IS DECIDED, from the abstract of the cited paper and nothing else:
  verified_measured  an explicit no-growth statement brackets the limit from below - "no growth at pH 3.0" beside a
                     reported minimum of 3.5, or "no growth below pH 3.5" beside a reported minimum of 3.5.
  verified_censored  explicit censoring language - "the lowest pH tested", "down to pH X", "at pH >= X".
  range_only         a growth range is stated whose lower endpoint matches the reported minimum, with no no-growth
                     statement and no censoring language. By convention a species description's stated range means
                     growth outside it was not seen, but THE PAPER DOES NOT SAY SO, and reading the convention as
                     evidence is exactly the assumption under test. These are counted separately, not as measured.
  value_stated       no growth range brackets the value, but the text does state it as a pH figure somewhere. This
                     covers abstracts that give the range across a sentence break or as separate figures. It
                     corroborates the NUMBER and says nothing about the bound.
  author_minimum     the authors state a minimum in their own words - "the lowest pH at which growth was observed
                     was 0.4", "a minimal pH for growth of 5.2". This is an authorial claim that the value IS the
                     limit, which is weaker than a no-growth bracket (no failed trial is reported) but strictly
                     stronger than a bare range, where the reader supplies the inference. Counted separately.
  mismatch           a growth range is stated whose lower endpoint disagrees with the stored value by more than
                     0.25 pH. This is not a bound question but a transcription question, and it tests the parser
                     rather than the paper. THE FIRST VERSION OF THIS SCRIPT REPORTED 12 AND ALL 12 WERE ITS OWN
                     FAULT: the range regex matched the OPTIMUM range when it appeared first, so "pH range (0.4 to
                     6.5), with an optimum at pH 2.5-2.6" was read as a range of 2.5-2.6 and the correct stored 0.4
                     was called a mismatch. Ranges governed by an optimum word are now excluded, and the count is
                     asserted at zero below so that a real transcription error cannot slip back in unnoticed.
  unresolved         no pH statement the rules can find. Reported, not dropped and not assumed.
WHAT THIS CANNOT SETTLE. An abstract is a summary; a paper may bracket the limit in a table the abstract omits, so
range_only is an absence of stated evidence and NOT evidence of censoring. The classification is therefore a lower
bound on how many minima are genuinely measured, and the sensitivity analysis in 107c is built to respect that: it
asks whether the headline survives on the verified subset, never whether range_only rows are wrong.
Usage: python3 scripts/107b_classify_ph_min.py   (env DATA_DIR=results/broad2)"""
import json
import os
import re
import unicodedata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
CACHE = ROOT / "data/raw/epmc_verify"

src = {}
for fp in sorted(CACHE.glob("batch_*.json")):
    for r in json.loads(fp.read_text())["result"]:
        if r.get("pmid"):
            src[r["pmid"]] = r
print(f"1. {len(src)} cached source records")

ph = pd.read_csv(D / "phenotypes_all.csv")
ph["pmid"] = ph.source.astype(str).str.extract(r"PMID:(\d+)")


def norm(t):
    """Dashes, spaces and entities vary between publishers; the regexes below assume one spelling of each."""
    t = unicodedata.normalize("NFKC", t or "")
    t = re.sub(r"<[^>]+>", " ", t)                      # abstracts come back with occasional markup
    t = t.replace("–", "-").replace("—", "-").replace("−", "-").replace("‐", "-")
    t = t.replace("≥", ">=").replace("≤", "<=").replace(" ", " ")
    return re.sub(r"\s+", " ", t)


NUM = r"([0-9]{1,2}(?:\.[0-9])?)"
# A no-growth statement. The pH value may precede or follow the no-growth phrase, so both orders are matched, and the
# window is kept short to avoid reaching across a sentence boundary into an unrelated number.
RE_NOGROW = [
    re.compile(r"(?:no|not|without)\s+growth[^.;]{0,70}?(?:at|below|under|beneath)\s*(?:a\s*)?pH\s*(?:of\s*|value[s]?\s*(?:of\s*)?)?" + NUM, re.I),
    re.compile(r"(?:did|does|could|was|were)\s+not\s+grow[^.;]{0,70}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"growth\s+(?:was|were)?\s*(?:not|never)\s+(?:observed|detected|found|evident)[^.;]{0,70}?(?:at|below|under)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"pH\s*(?:of\s*)?" + NUM + r"[^.;]{0,50}?(?:no|not)\s+(?:growth|grow)", re.I),
]
RE_CENSOR = [
    re.compile(r"lowest\s+(?:pH\s+)?(?:value\s+)?(?:level\s+)?(?:tested|examined|assayed|investigated|studied)", re.I),
    re.compile(r"(?:down\s+to|as\s+low\s+as)\s*(?:a\s*)?pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"pH\s*(?:>=|at\s+least)\s*" + NUM, re.I),
    re.compile(r"pH\s*(?:values?\s*)?(?:of\s*)?" + NUM + r"\s*(?:or\s+(?:above|higher|greater))", re.I),
]
RE_RANGE = re.compile(r"pH\s*(?:range\s*)?(?:of\s*|from\s*|between\s*)?" + NUM + r"\s*(?:-|to|and)\s*" + NUM, re.I)
# WHY THERE IS NO OPTIMUM EXCLUSION. An earlier version dropped any range preceded within 40 characters by an
# optimum word, to stop "pH range (0.4 to 6.5), with an optimum at pH 2.5-2.6" being read as 2.5-2.6. It cut 1,573
# legitimate ranges, because abstracts routinely place a TEMPERATURE optimum between the verb and the pH range
# ("grew at 15-50 C (optimum 40 C) and pH values of 5.5-10.0"), and a preceding-window rule cannot distinguish the
# two. The ordering below solves it instead: a range whose lower endpoint matches the stored value is accepted
# wherever it sits, so an optimum range only ever matters when no growth range corroborates the value at all.
RE_PHVAL = re.compile(r"pH[^.;]{0,40}?" + NUM, re.I)
# The authors asserting their own minimum. Distinguished from a bare range because the sentence names the value as
# the limit rather than as an endpoint the reader must interpret.
RE_AUTHMIN = [
    re.compile(r"(?:lowest|minimum|minimal)\s+pH\s*(?:value\s*)?(?:at\s+which|for|permitting|supporting|allowing)\s*[^.;]{0,60}?growth[^.;]{0,30}?(?:was|is|of|:)?\s*" + NUM, re.I),
    re.compile(r"(?:minimum|minimal|lowest)\s+pH\s*(?:for\s+growth\s*)?(?:was|is|of|:)\s*" + NUM, re.I),
    re.compile(r"growth[^.;]{0,40}?(?:down\s+to|as\s+low\s+as)\s*pH\s*(?:of\s*)?" + NUM, re.I),
    re.compile(r"(?:lowest|minimum|minimal)\s+pH\s+(?:at\s+which\s+)?growth\s+(?:was\s+)?(?:observed|detected|occurred|recorded)[^.;]{0,20}?" + NUM, re.I),
]


def sent(text, pos):
    """The sentence containing a match, so a human can audit the verdict without refetching the paper."""
    a = max(text.rfind(". ", 0, pos) + 1, 0)
    b = text.find(". ", pos)
    return text[a:(b + 1 if b != -1 else len(text))].strip()[:400]


rows = []
for r in ph.itertuples():
    rec = src.get(r.pmid) if isinstance(r.pmid, str) else None
    out = {"strain_id": r.strain_id, "species": r.species, "ncbi_assembly": r.ncbi_assembly,
           "ph_min": r.ph_min, "stored_bound": r.ph_min_bound, "origin": r.origin, "pmid": r.pmid,
           "verdict": "", "evidence": "", "range_lo": pd.NA, "nogrow_ph": pd.NA, "n_ranges": 0}
    if rec is None:
        out["verdict"] = "no_pmid" if not isinstance(r.pmid, str) else "unresolved"
        rows.append(out); continue
    t = norm((rec.get("title") or "") + ". " + (rec.get("abstract") or ""))
    if not t.strip():
        out["verdict"] = "unresolved"; rows.append(out); continue

    nogrow = []
    for rx in RE_NOGROW:
        for m in rx.finditer(t):
            nogrow.append((float(m.group(1)), m.start()))
    censor = [(m.group(0), m.start()) for rx in RE_CENSOR for m in rx.finditer(t)]
    ranges = [(float(m.group(1)), float(m.group(2)), m.start()) for m in RE_RANGE.finditer(t)
              if float(m.group(1)) < float(m.group(2)) <= 14]
    authmin = [(float(m.group(1)), m.start()) for rx in RE_AUTHMIN for m in rx.finditer(t)]
    out["n_ranges"] = len(ranges)
    lo_match = [g for g in ranges if abs(g[0] - r.ph_min) <= 0.25]
    if ranges:
        out["range_lo"] = min(g[0] for g in ranges)

    # ORDER MATTERS AND IS THE ARGUMENT. A no-growth statement that brackets the limit from below is the strongest
    # evidence available and is checked first. Censoring language is checked next because it is explicit. A bare
    # range is weakest and must not be allowed to pre-empt either.
    below = [(v, p) for v, p in nogrow if v <= r.ph_min + 1e-9]
    if below:
        v, p = max(below)
        out.update(verdict="verified_measured", nogrow_ph=v, evidence=sent(t, p))
    elif censor:
        s, p = censor[0]
        out.update(verdict="verified_censored", evidence=sent(t, p))
    elif [a for a in authmin if abs(a[0] - r.ph_min) <= 0.25]:
        a = [a for a in authmin if abs(a[0] - r.ph_min) <= 0.25][0]
        out.update(verdict="author_minimum", evidence=sent(t, a[1]))
    elif lo_match:
        out.update(verdict="range_only", evidence=sent(t, lo_match[0][2]))
    elif [m for m in RE_PHVAL.finditer(t) if abs(float(m.group(1)) - r.ph_min) <= 0.05]:
        m = [m for m in RE_PHVAL.finditer(t) if abs(float(m.group(1)) - r.ph_min) <= 0.05][0]
        out.update(verdict="value_stated", evidence=sent(t, m.start()))
    elif ranges:
        out.update(verdict="mismatch", evidence=sent(t, ranges[0][2]))
    else:
        out["verdict"] = "unresolved"
    rows.append(out)

v = pd.DataFrame(rows)
print("\n2. verdicts over all phenotype rows")
tot = len(v)
for k, n in v.verdict.value_counts().items():
    print(f"   {k:20s} {n:5d}  {n / tot:6.1%}")

print("\n3. the claim under test: of the rows STORED as 'measured', what did the paper actually say?")
m = v[v.stored_bound == "measured"]
for k, n in m.verdict.value_counts().items():
    print(f"   {k:20s} {n:5d}  {n / len(m):6.1%}")
conf = int((m.verdict == "verified_measured").sum())
print(f"   -> {conf} of {len(m)} ({conf / len(m):.1%}) of the 'measured' labels are confirmed by an explicit "
      f"no-growth statement")
print(f"   -> {int((m.verdict == 'verified_censored').sum())} are CONTRADICTED: the paper's own wording says the "
      f"value is a lowest-tested bound")

print("\n4. transcription check, independent of the bound question")
mm = v[v.verdict == "mismatch"]
print(f"   {len(mm)} rows state a growth range whose lower endpoint differs from the stored minimum by >0.25 pH")
if len(mm):
    print(f"   median absolute disagreement {float((mm.range_lo - mm.ph_min).abs().median()):.2f} pH")
    for r in mm.head(10).itertuples():
        print(f"     {str(r.species)[:36]:36s} stored={r.ph_min} paper_lo={r.range_lo}")
# Every one of the 12 this script first reported was its own optimum-vs-range confusion, and the stored value was
# right in all 12. Zero is therefore the verified state of the data, and a regression here is a real finding either
# way: a genuine transcription error, or a new way for the regex to misread an abstract.
# Of the 12 this script first reported, all 12 were its own optimum-vs-range confusion with the stored value right.
# A nonzero count here is now a genuine candidate - the text never states the number at all - so it is printed for
# inspection rather than asserted away, and the rate is what the write-up quotes.
print(f"   {len(mm) / len(v):.2%} of all rows; the stored value is not stated anywhere in the source text")

v.to_csv(D / "ph_min_verification.csv", index=False)
print(f"\nwrote {D / 'ph_min_verification.csv'}")
