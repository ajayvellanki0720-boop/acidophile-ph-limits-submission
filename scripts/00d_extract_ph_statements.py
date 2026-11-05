"""Extract pH/temperature statements from cached Europe PMC abstracts into a REVIEW table
(data/curated/epmc_ph_review.csv). Regex output is a first pass: a person (or Claude, reading each row)
confirms taxon, strain, values, and bound type before anything enters the phenotype table."""
import glob
import html
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N = r"(\d+(?:\.\d+)?)"
SEP = r"\s*(?:–|-|—|−|to|and)\s*"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth))?\s*" + N + SEP + N)
OPT = [re.compile(r"optim\w*[^.;]{0,40}?pH\s*" + N + r"(?:" + SEP + N + r")?"),
       re.compile(r"pH\s*optim\w*\s*(?:of|at|:|,|was|were)?\s*" + N + r"(?:" + SEP + N + r")?")]
TEMP = re.compile(N + SEP + N + r"\s*°\s*C(?:[^.]{0,20}?optim\w*[,:]?\s*(?:at\s*)?" + N + r")?")
HEDGE = re.compile(r"at least|not tested|no lower|lowest|as low as|down to|as high as|up to|above|below|not determined|only", re.I)
SPECIES = re.compile(r"([A-Z][a-z]{3,}) ([a-z]{3,})\s*(?:subsp\. [a-z]+\s*)?(?:sp\. nov|comb\. nov|nov\.)")
TITLE_BINOMIAL = re.compile(r"^(?:Description of |Isolation and characterization of |Emended description of )?(?:Candidatus )?([A-Z][a-z]{3,}) ([a-z]{4,})")
GENUS = re.compile(r"([A-Z][a-z]{3,})\s+gen\. nov")

rows = []
seen = set()
for f in sorted(glob.glob(str(ROOT / "data/raw/europepmc/*.json"))):
    for r in json.load(open(f)):
        key = r["id"]
        if key in seen or not r.get("abstractText"):
            continue
        seen.add(key)
        title = html.unescape(re.sub(r"<[^>]+>", "", r.get("title") or ""))
        ab = html.unescape(re.sub(r"<[^>]+>", "", r["abstractText"]))
        m = TITLE_BINOMIAL.search(" ".join(title.split())) or SPECIES.search(title) or SPECIES.search(ab)
        g = GENUS.search(title)
        sents = [s for s in re.split(r"(?<=[a-z\)\]])\.\s+(?=[A-Z])", ab) if "pH" in s]
        rng = RANGE.search(" ".join(sents))
        if not rng:
            continue
        txt = " ".join(sents)
        opt = next((o.search(txt) for o in OPT if o.search(txt)), None)
        tmp = TEMP.search(ab)
        rows.append({
            "pmid": r.get("pmid") or key, "doi": r.get("doi"), "year": r.get("pubYear"), "journal": r.get("journalTitle"),
            "taxon": f"{m.group(1)} {m.group(2)}" if m else (f"{g.group(1)} (gen. nov.)" if g else None),
            "title": title[:140], "ph_lo": float(rng.group(1)), "ph_hi": float(rng.group(2)),
            "ph_opt_lo": float(opt.group(1)) if opt else None,
            "ph_opt_hi": float(opt.group(2)) if opt and opt.group(2) else None,
            "temp_lo": float(tmp.group(1)) if tmp else None, "temp_hi": float(tmp.group(2)) if tmp else None,
            "temp_opt": float(tmp.group(3)) if tmp and tmp.group(3) else None,
            "hedge": bool(HEDGE.search(txt)), "statement": txt[:420],
        })
out = pd.DataFrame(rows).drop_duplicates("pmid")
out = out[out.ph_lo <= out.ph_hi]
out.to_csv(ROOT / "data/curated/epmc_ph_review.csv", index=False)
print(len(out), "abstracts with a pH range;", out.taxon.notna().sum(), "with a taxon;", out.ph_opt_lo.notna().sum(), "with optimum;", out.hedge.sum(), "hedged")
