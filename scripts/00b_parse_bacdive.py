"""Parse raw BacDive JSON into a CANDIDATE phenotype table (data/curated/bacdive_candidates.csv).

BacDive reports pH ranges without saying whether growth failed beyond them, so every ph_min/ph_max
here gets bound type `unknown` and confidence `low` until a primary source is checked. Nothing in this
file can enter the primary model set without that verification. Strains with no pH are dropped."""
import glob
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.phenotype import REQUIRED_COLUMNS

ROOT = Path(__file__).resolve().parents[1]


def as_list(x):
    return x if isinstance(x, list) else ([x] if x else [])


def nums(s):
    return [float(t) for t in re.findall(r"\d+(?:\.\d+)?", str(s))]


rows = []
for f in sorted(glob.glob(str(ROOT / "data/raw/bacdive/*.json"))):
    for bid, v in json.load(open(f)).items():
        c = v.get("Culture and growth conditions", {})
        ph = as_list(c.get("culture pH"))
        if not ph:
            continue
        lo = hi = opt = None
        for e in ph:
            n = nums(e.get("pH"))
            if not n or e.get("ability") == "negative":
                continue
            t = e.get("type")
            if t == "optimum":
                opt = sum(n) / len(n)
            elif t == "growth" and len(set(n)) > 1:  # a single value is a growth point, not a range
                lo, hi = min(n), max(n)
            elif t == "minimum":
                lo = min(n)
            elif t == "maximum":
                hi = max(n)
        taxo = v["Name and taxonomic classification"]
        temps = as_list(c.get("culture temp"))
        topt = [nums(e.get("temperature"))[0] for e in temps if e.get("type") == "optimum" and nums(e.get("temperature"))]
        tgr = [x for e in temps if e.get("type") == "growth" for x in nums(e.get("temperature"))]
        acc = [g["accession"] for g in as_list(v.get("Sequence information", {}).get("Genome sequences"))
               if g.get("database") == "ncbi" and str(g.get("accession", "")).startswith(("GCA_", "GCF_"))]
        rows.append({
            "strain_id": f"BD{bid}", "bacdive_id": bid, "ncbi_assembly": acc[0] if acc else None,
            "n_assemblies": len(acc), "domain": taxo.get("domain"), "phylum": taxo.get("phylum"),
            "class": taxo.get("class"), "family": taxo.get("family"), "genus": taxo.get("genus"),
            "species": taxo.get("species"), "strain_designation": taxo.get("strain designation"),
            "ph_min": lo, "ph_min_bound": "unknown", "ph_opt": opt, "ph_max": hi, "ph_max_bound": "unknown",
            "acid_identity": "unspecified", "organic_acid": None, "organic_acid_mM": None, "ph_measured_at_c": None,
            "temp_opt_c": topt[0] if topt else None, "temp_min_c": min(tgr) if tgr else None,
            "temp_max_c": max(tgr) if tgr else None, "oxygen": None, "medium": None, "metabolism": None,
            "metal_note": None, "source": f"BacDive {bid}", "confidence": "low",
        })
out = pd.DataFrame(rows)
for col in REQUIRED_COLUMNS:
    if col not in out:
        out[col] = None
out.to_csv(ROOT / "data/curated/bacdive_candidates.csv", index=False)
print(f"{len(out)} strains with pH | with min {out.ph_min.notna().sum()} | with opt {out.ph_opt.notna().sum()} | "
      f"with assembly {out.ncbi_assembly.notna().sum()} | min+opt+assembly {(out.ph_min.notna() & out.ph_opt.notna() & out.ncbi_assembly.notna()).sum()}")
print(out.groupby("genus").size().sort_values(ascending=False).to_string())
