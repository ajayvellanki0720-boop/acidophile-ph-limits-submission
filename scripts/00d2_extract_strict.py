"""Stricter second-pass extraction from all cached Europe PMC abstracts -> data/curated/epmc_review_v2.csv.
Keeps only species/genus descriptions with an acid-relevant growth pH range and a single named taxon, and marks
`clean` when every field parsed unambiguously. Rows not clean still need manual reading. Already-reviewed PMIDs
(in 00f CURATED/EXCLUDED) are skipped."""
import glob
import html
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N = r"(\d+(?:\.\d+)?)"
SEP = r"\s*(?:–|-|—|−|to|and)\s*"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + SEP + N)
BAD_CTX = re.compile(r"isolat|sample|collected|habitat|spring|drainage|soil|sediment|environment|originally|measured|in situ|bog|lake|water of", re.I)
OPT_AFTER = re.compile(r"[\(,;]\s*(?:pH\s*)?optim\w*[,:]?\s*(?:growth\s*)?(?:at|of|around|about|~|was)?\s*(?:pH\s*)?" + N + r"(?:" + SEP + N + r")?")
OPT_PHRASE = re.compile(r"(?:optim\w*\s*(?:growth\s*)?(?:pH|at pH|was pH|pH of)\s*" + N + r"(?:" + SEP + N + r")?)|(?:pH\s*optim\w*\s*(?:of|at|:|was)?\s*" + N + r"(?:" + SEP + N + r")?)", re.I)
DESC = re.compile(r"sp\. nov|gen\. nov|comb\. nov|emended|description of|reclassif|novel (species|genus)|proposal", re.I)
HEDGE = re.compile(r"at least|not tested|lowest|as low as|down to|as high as|up to|only|also (grew|occurred)|good growth|initial pH", re.I)

src = Path(ROOT / "scripts/00f_curate_literature.py").read_text()
done = set(re.findall(r'^\s+"([\w]+)": [(\"]', src, re.M))
rows, seen = [], set(done)
for f in sorted(glob.glob(str(ROOT / "data/raw/europepmc/*.json"))):
    for r in json.load(open(f)):
        pid = r.get("pmid") or r["id"]
        if pid in seen or not r.get("abstractText"):
            continue
        seen.add(pid)
        title = " ".join(html.unescape(re.sub(r"<[^>]+>", "", r.get("title") or "")).split())
        ab = " ".join(html.unescape(re.sub(r"<[^>]+>", "", r["abstractText"])).split())
        if not DESC.search(title + " " + ab[:300]):
            continue
        sents = [s for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", ab) if "pH" in s]
        cands = []
        for s in sents:
            for m in RANGE.finditer(s):
                lo, hi = float(m.group(1)), float(m.group(2))
                if lo < hi and not BAD_CTX.search(s[max(0, m.start() - 45):m.start()]) and re.search(r"grow|growth|grew", s, re.I):
                    cands.append((hi - lo, lo, hi, s, m))
        if not cands:
            continue
        _, lo, hi, s, m = max(cands, key=lambda c: c[0])
        if lo > 5.5 and hi < 7:       # clearly neutral-only
            continue
        tail = s[m.end():m.end() + 70]
        o = OPT_AFTER.match(tail.lstrip()) if False else OPT_AFTER.search(tail)
        if not o:
            o = OPT_PHRASE.search(s)
        opt = None
        if o:
            span = s[o.start():o.end()]
            if not re.search(r"°|degrees|NaCl|temperature", span):
                v = [float(g) for g in o.groups() if g]
                opt = sum(v) / len(v) if v else None
        tm = re.match(r"^(?:Description of |Isolation and characterization of |Emended description of )?(?:Candidatus )?([A-Z][a-z]{3,}) ([a-z]{4,})", title)
        n_taxa = len(re.findall(r"sp\. nov|comb\. nov", title))
        multi = n_taxa > 1 or bool(re.search(r"\bsp\. nov\.?,? [A-Z][a-z]+ [a-z]+ sp\. nov|and [A-Z][a-z]+ [a-z]+ sp\. nov", title))
        hedged = bool(HEDGE.search(s))
        clean = bool(tm) and not multi and not hedged and opt is not None and lo <= 5.5 and (hi - lo) >= 1.0
        rows.append({"pmid": pid, "year": r.get("pubYear"), "taxon": f"{tm.group(1)} {tm.group(2)}" if tm else None,
                     "multi_taxon_title": multi, "ph_lo": lo, "ph_hi": hi, "ph_opt": opt, "hedged": hedged, "clean": clean,
                     "title": title[:130], "statement": s[max(0, m.start() - 60):m.start() + 240]})
out = pd.DataFrame(rows)
out.to_csv(ROOT / "data/curated/epmc_review_v2.csv", index=False)
print(len(out), "rows |", out.clean.sum(), "clean |", (~out.clean & out.taxon.notna() & ~out.multi_taxon_title).sum(), "need manual |",
      out.multi_taxon_title.sum(), "multi-taxon")
