"""Rule-based extraction of growth-pH range, optimum and taxon from the broad species-description abstracts, with VALIDATION against the
hand-curated literature rows (data/curated/literature_curated.csv). Output: data/curated/broad_candidates.csv (all parsed rows, flagged) and
results/full/parser_validation.csv. Only rows flagged `clean` are eligible for automatic use; everything else needs reading.
Usage: python3 scripts/23_parse_broad.py"""
import glob
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N = r"(\d+(?:\.\d+)?)"
SEP = r"\s*(?:–|-|—|−|to|and)\s*"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + SEP + N)
BAD_CTX = re.compile(r"isolat|sample|collected|habitat|spring|drainage|soil|sediment|environment|originally|measured|in situ|bog|lake|water of|enzyme|activity|stable|stability|protease|lipase|amylase|cellulase", re.I)
OPT_AFTER = re.compile(r"[\(,;]\s*(?:pH\s*)?optim\w*[,:]?\s*(?:growth\s*)?(?:at|of|around|about|~|was)?\s*(?:pH\s*)?" + N + r"(?:" + SEP + N + r")?")
OPT_PHRASE = re.compile(r"(?:optim\w*\s*(?:growth\s*)?(?:pH|at pH|was pH|pH of)\s*" + N + r"(?:" + SEP + N + r")?)|(?:pH\s*optim\w*\s*(?:of|at|:|was)?\s*" + N + r"(?:" + SEP + N + r")?)", re.I)
HEDGE = re.compile(r"at least|not tested|lowest|as low as|down to|as high as|up to|only|also (grew|occurred)|good growth|initial pH|tolerat", re.I)
MULTI = re.compile(r"respectively|the former|the latter|whereas|while strain|strains? [A-Z0-9-]+T? (?:and|,) [A-Z0-9-]+T? (?:grew|were able)", re.I)
TITLE_BIN = re.compile(r"^(?:Description of |Isolation and characterization of |Emended description of |Proposal of )?(?:Candidatus |“Candidatus |\"Candidatus )?([A-Z][a-z]{3,}) ([a-z]{4,})")


def parse(rec):
    title, ab = rec["title"], rec["abstract"]
    sents = [s for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", ab) if "pH" in s]
    cands = []
    for s in sents:
        for m in RANGE.finditer(s):
            lo, hi = float(m.group(1)), float(m.group(2))
            if lo < hi and hi <= 14 and not BAD_CTX.search(s[max(0, m.start() - 45):m.start()]) and re.search(r"grow|growth|grew|optimum|range", s, re.I):
                cands.append((hi - lo, lo, hi, s, m))
    if not cands:
        return None
    _, lo, hi, s, m = max(cands, key=lambda c: c[0])
    o = OPT_AFTER.search(s[m.end():m.end() + 70]) or OPT_PHRASE.search(s)
    opt = None
    if o:
        span = s[o.start():o.end()]
        after = s[o.end() - 0:o.end() + 8]            # a number followed by a unit is a temperature / salinity, not a pH
        if not re.search(r"°|degrees|NaCl|temperature", span) and not re.match(r"\s*(?:°|º|degrees|%|\s?C\b|mM|g\b|M\b)", after):
            v = [float(g) for g in o.groups() if g]; opt = sum(v) / len(v) if v else None
            if opt is not None and not (lo - 1e-9 <= opt <= hi + 1e-9):     # the optimum must lie inside the growth range
                opt = None
    tm = TITLE_BIN.search(title)
    n_taxa = len(re.findall(r"sp\. nov|comb\. nov", title)); multi_title = n_taxa > 1 or bool(re.search(r" and [A-Z][a-z]+ [a-z]+ sp\. nov", title))
    multi_text = bool(MULTI.search(s)); hedged = bool(HEDGE.search(s))
    clean = bool(tm) and not multi_title and not multi_text and not hedged and opt is not None and (hi - lo) >= 1.0 and lo <= hi and 0 <= lo <= 12
    usable_noopt = bool(tm) and not multi_title and not multi_text and not hedged and (hi - lo) >= 1.0 and 0 <= lo <= 12     # as clean, but no optimum required
    return {"pmid": rec.get("pmid") or rec["id"], "doi": rec.get("doi"), "year": rec["year"], "taxon": f"{tm.group(1)} {tm.group(2)}" if tm else None, "ph_min": lo, "ph_opt": opt, "ph_max": hi,
            "multi_title": multi_title, "multi_text": multi_text, "hedged": hedged, "clean": clean, "usable_noopt": usable_noopt, "title": title[:130], "statement": s[max(0, m.start() - 60):m.start() + 220]}


rows, seen = [], set()
for f in sorted(glob.glob(str(ROOT / "data/raw/europepmc_broad/*.json"))):
    for rec in json.load(open(f)):
        if rec["id"] in seen:
            continue
        seen.add(rec["id"]); r = parse(rec)
        if r:
            rows.append(r)
B = pd.DataFrame(rows); B.to_csv(ROOT / "data/curated/broad_candidates.csv", index=False)
print(len(B), "parsed rows |", int(B.clean.sum()), "clean |", int(B.multi_title.sum()), "multi-taxon titles |", int(B.hedged.sum()), "hedged")
# validation against the hand-curated rows (matched by PMID)
cur = pd.read_csv(ROOT / "data/curated/literature_curated.csv", dtype={"pmid": str})
cur = cur[cur.bound == "measured"].drop_duplicates("pmid", keep="first")
V = B.astype({"pmid": str}).merge(cur[["pmid", "taxon", "ph_min", "ph_opt", "ph_max"]], on="pmid", suffixes=("_auto", "_hand"))
V["min_ok"] = (V.ph_min_auto - V.ph_min_hand).abs() <= 0.05; V["max_ok"] = (V.ph_max_auto - V.ph_max_hand).abs() <= 0.05
V["opt_ok"] = ((V.ph_opt_auto - V.ph_opt_hand).abs() <= 0.3) | (V.ph_opt_auto.isna() & V.ph_opt_hand.isna())
V["taxon_ok"] = V.taxon_auto == V.taxon_hand
V["all_ok"] = V.min_ok & V.max_ok & V.opt_ok & V.taxon_ok
V.to_csv(ROOT / "results/full/parser_validation.csv", index=False)
for lab, sub in (("all matched", V), ("flagged clean", V[V.clean])):
    if len(sub):
        print(f"validation, {lab}: n={len(sub)} | pH_min correct {sub.min_ok.mean():.2f} | pH_max {sub.max_ok.mean():.2f} | optimum {sub.opt_ok.mean():.2f} | taxon {sub.taxon_ok.mean():.2f} | all four {sub.all_ok.mean():.2f}")
