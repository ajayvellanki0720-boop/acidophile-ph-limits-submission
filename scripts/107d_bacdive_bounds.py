"""Harvest explicit growth / no-growth pH records from BacDive for every genus in the panel -> DATA_DIR/ph_min_bacdive.csv
WHY THIS IS THE ONE ROUTE LEFT. 107b established that the source abstracts cannot settle whether a minimum is a
measured limit or a lowest-tested bound: 94.9% state a pH range, but only 1.1% contain "no growth" at all and none
contains "lowest pH tested", because a species description states the range and leaves the bracketing to its tables.
107c established that those tables are reachable for 2.5% of sources, open access being the limit, and that route is
exhausted. BacDive is different in kind. It does not store prose; it stores a typed record per observation:
    {"ability": "no",       "type": "growth",  "pH": "2.5"}   <- an explicit FAILED trial at pH 2.5
    {"ability": "positive", "type": "growth",  "pH": "1-3"}   <- a growth range
    {"ability": "positive", "type": "minimum", "pH": "1.5"}   <- a curator-assigned minimum
An `ability: no` record for a pH at or below the stored minimum is exactly the bracketing evidence the abstracts
withhold, and it is structured rather than inferred. In the 31 genera already cached, 8 of 79 culture-pH entries
carry `ability: no`, so the route is worth running across all 1,366 genera of the panel.
WHAT THE PARTIAL RUN ALREADY SHOWS, recorded here because it decides how the result must be read. Of the first
21 minima BacDive brackets, every one is a neutrophile: stored minima run 4.1 to 7.5 and cluster at 6-7. Not one
acidophile is bracketed. The reason is not curation quality but what experimenters do - a failed trial at pH 4 or 5
gets recorded for an organism nobody expected to grow there, while for an organism whose minimum is pH 0.4 nobody
records a failure below it. So this route supplies bound evidence precisely where censoring was never a threat
(organisms with pH units of headroom) and supplies none where it would bite (the acid-tolerant tail the study is
about). That asymmetry must be stated with any count from this script, or the count will read as reassurance it
cannot give.
WHAT IT IS NOT. BacDive is a curation of the same primary literature, so it is a better-structured reading of the
same papers, not an independent measurement. Agreement with the stored value is therefore weak corroboration; a
no-growth record is strong, because it reports a trial that the abstract did not mention. Where BacDive and the
stored value disagree the row is flagged, not overwritten: scripts/00g already prefers literature over BacDive by
design, and this script exists to audit that preference rather than to invert it.
Usage: python3 scripts/107d_bacdive_bounds.py   (env DATA_DIR=results/broad2)"""
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
RAW = ROOT / "data/raw/bacdive_panel"
RAW.mkdir(parents=True, exist_ok=True)
UA = "acidophile-ph-limits research (open API; polite delay)"
DELAY = 0.6


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                return None                      # genus simply absent from BacDive; a fact, not an error
            if i == tries - 1:
                return {"__error__": f"HTTP {e.code}"}
            time.sleep(3 * (i + 1))
        except Exception:
            if i == tries - 1:
                return {"__error__": "network"}
            time.sleep(3 * (i + 1))


ph = pd.read_csv(D / "phenotypes_all.csv")
genera = sorted(ph.species.astype(str).str.split().str[0].unique())
print(f"1. {len(genera)} genera across {len(ph)} panel rows", flush=True)

done = {p.stem for p in RAW.glob("*.json")}
todo = [g for g in genera if g not in done]
print(f"2. {len(done)} already cached, {len(todo)} to fetch "
      f"(~{len(todo) * DELAY * 2.2 / 60:.0f} min at {DELAY}s/call)", flush=True)

for i, g in enumerate(todo, 1):
    ids_page = get(f"https://api.bacdive.dsmz.de/taxon/{urllib.request.quote(g)}")
    time.sleep(DELAY)
    if not ids_page or "__error__" in ids_page:
        (RAW / f"{g}.json").write_text(json.dumps({"genus": g, "records": {},
                                                   "note": (ids_page or {}).get("__error__", "absent")}))
        continue
    ids = list(ids_page.get("results") or [])
    recs = {}
    for c in range(0, len(ids), 100):
        chunk = ids[c:c + 100]
        r = get("https://api.bacdive.dsmz.de/fetch/" + ";".join(map(str, chunk)))
        time.sleep(DELAY)
        if r and "__error__" not in r:
            recs.update(r.get("results") or {})
    (RAW / f"{g}.json").write_text(json.dumps({"genus": g, "records": recs}))
    if i % 50 == 0:
        print(f"   {i}/{len(todo)} genera, last={g}, {len(recs)} strains", flush=True)

print("3. extracting typed pH observations", flush=True)


def vals(s):
    """A BacDive pH field is a scalar or a range, written with any dash. Returns every number it contains."""
    return [float(x) for x in re.findall(r"[0-9]{1,2}(?:\.[0-9]+)?", str(s).replace("–", "-"))]


obs = []
for fp in sorted(RAW.glob("*.json")):
    blob = json.loads(fp.read_text())
    for bid, rec in (blob.get("records") or {}).items():
        if not isinstance(rec, dict):
            continue
        name = (((rec.get("Name and taxonomic classification") or {}).get("species"))
                or ((rec.get("Name and taxonomic classification") or {}).get("full scientific name")) or "")
        name = re.sub(r"<[^>]+>", "", str(name)).strip()
        cp = ((rec.get("Culture and growth conditions") or {}).get("culture pH"))
        if not cp:
            continue
        for e in (cp if isinstance(cp, list) else [cp]):
            if not isinstance(e, dict):
                continue
            v = vals(e.get("pH"))
            if not v:
                continue
            obs.append({"bacdive_id": bid, "bacdive_species": name, "ability": e.get("ability"),
                        "type": e.get("type"), "ph_lo": min(v), "ph_hi": max(v), "ref": e.get("@ref")})
o = pd.DataFrame(obs)
print(f"   {len(o)} typed pH observations over {o.bacdive_species.nunique() if len(o) else 0} BacDive species",
      flush=True)
if len(o):
    print("   ability x type:\n" + pd.crosstab(o.ability, o.type).to_string())
o.to_csv(D / "ph_min_bacdive_observations.csv", index=False)

print("\n4. matching to the panel and deciding the bound", flush=True)
rows = []
by_sp = {k: g for k, g in o.groupby("bacdive_species")} if len(o) else {}
for r in ph.itertuples():
    g = by_sp.get(str(r.species))
    out = {"species": r.species, "ph_min": r.ph_min, "stored_bound": r.ph_min_bound,
           "bacdive_verdict": "no_bacdive_record", "bacdive_evidence": "", "n_obs": 0}
    if g is not None and len(g):
        out["n_obs"] = len(g)
        no = g[(g.ability.astype(str) == "no") & (g.type.astype(str) == "growth")]
        brack = no[no.ph_lo < r.ph_min - 1e-9]        # strictly below: a failed trial beneath the stored limit
        contra = no[(no.ph_lo - r.ph_min).abs() <= 1e-9]  # at the limit itself: the record contradicts the value
        gro = g[(g.ability.astype(str) == "positive") & (g.type.astype(str) == "growth")]
        mini = g[g.type.astype(str) == "minimum"]
        if len(brack):
            out.update(bacdive_verdict="verified_measured",
                       bacdive_evidence=f"no growth recorded at pH {brack.ph_lo.max():g}, below the stored minimum")
        elif len(contra):
            out.update(bacdive_verdict="contradicts_stored_minimum",
                       bacdive_evidence=f"no growth recorded AT pH {contra.ph_lo.iloc[0]:g}, the stored minimum itself")
        elif len(mini) and (mini.ph_lo - r.ph_min).abs().min() <= 0.25:
            out.update(bacdive_verdict="curator_minimum",
                       bacdive_evidence=f"BacDive minimum pH {mini.ph_lo.iloc[0]:g}")
        elif len(gro) and (gro.ph_lo - r.ph_min).abs().min() <= 0.25:
            out.update(bacdive_verdict="range_only",
                       bacdive_evidence=f"growth range from pH {gro.ph_lo.min():g}")
        elif len(gro):
            out.update(bacdive_verdict="range_disagrees",
                       bacdive_evidence=f"BacDive lowest growth pH {gro.ph_lo.min():g} vs stored {r.ph_min:g}")
        else:
            out.update(bacdive_verdict="record_without_growth_ph")
    rows.append(out)
b = pd.DataFrame(rows)
for k, n in b.bacdive_verdict.value_counts().items():
    print(f"   {k:26s} {n:5d}  {n / len(b):6.1%}", flush=True)
b.to_csv(D / "ph_min_bacdive.csv", index=False)
print(f"\nwrote {D / 'ph_min_bacdive.csv'}", flush=True)
