"""Broader Europe PMC sweep: acid-related species descriptions by publication year (raw -> data/raw/europepmc/year_<y>.json)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

m = import_module("00c_fetch_europepmc")
ACID = ('(ABSTRACT:"acidophil*" OR ABSTRACT:"acidotoler*" OR ABSTRACT:"acidiphil*" OR ABSTRACT:"acidic" OR ABSTRACT:"acid mine" '
        'OR ABSTRACT:"low pH" OR ABSTRACT:"thermoacidophil*" OR ABSTRACT:"acid soil" OR ABSTRACT:"peat")')
for y in range(1990, 2027):
    f = m.RAW / f"year_{y}.json"
    if f.exists():
        continue
    q = f'(ABSTRACT:"sp. nov." OR ABSTRACT:"gen. nov." OR ABSTRACT:"comb. nov.") AND ABSTRACT:"pH" AND {ACID} AND PUB_YEAR:{y}'
    try:
        r = m.fetch(q)
    except Exception as e:
        print(y, "FAILED", e); continue
    f.write_text(json.dumps(r)); print(y, len(r), flush=True)
