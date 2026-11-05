"""Select mutagenesis targets from the per-position scan (37a) -> DATA_DIR/mutagenesis_targets.csv. Re-scores the saved table, so it needs no realignment.
The scan showed that exposed periplasmic positions are mostly NOT polymorphic between acidic and basic (only 4 of 250 carry each state in at least 15 % of bacteria): most
genomes have neither D/E nor K/R there, and the association is an aggregate over many positions rather than a swap at a few signature sites. Targets are therefore chosen to
change the aggregate safely rather than to copy a natural polymorphism:
  - the Escherichia coli reference carries D or E, so a reversal to K is possible;
  - at least 2 % of bacteria carry a basic residue at that position, so a basic state is known to be tolerated somewhere;
  - the position is not a documented loss-of-function or activity-reducing site (UniProt Mutagenesis features) and, for YidC, lies in 24-264, the region reported to be
    removable without lethality, so a charge change there tests electrostatics rather than folding;
  - positions are then ranked by their own association (acidic state vs pH_min, genome-wide composition controlled), with Benjamini-Hochberg q over the 165 positions where
    that association could be estimated. Individually significant positions are NOT required: the design changes the aggregate, and per-position ranking only orders the list.
Usage: python3 scripts/37b_select_targets.py   (env DATA_DIR=results/broad2)"""
import json
import os
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
ACC = {"YidC": "P25714", "NuoL": "P33607", "NuoN": "P0AFF0"}
SAFE_REGION = {"YidC": (24, 264)}   # UniProt: "can be removed without causing lethality"
R = pd.read_csv(D / "candidate_positions.csv")
bad = {}
for ref, acc in ACC.items():
    j = json.loads((ROOT / "data/raw/structures" / f"{acc}.json").read_text())
    s = set()
    for f in j["features"]:
        if f["type"] == "Mutagenesis":
            s.update(range(f["location"]["start"]["value"], f["location"]["end"]["value"] + 1))
    bad[ref] = s
R["q"] = pd.NA
msk = R.beta.notna()
R.loc[msk, "q"] = stats.false_discovery_control(R.loc[msk, "p"].to_numpy(), method="bh")
R["known_site"] = [p in bad.get(r, set()) for r, p in zip(R.ref, R.pos)]
R["in_safe_region"] = [not (r in SAFE_REGION) or (SAFE_REGION[r][0] <= p <= SAFE_REGION[r][1]) for r, p in zip(R.ref, R.pos)]
R["basic_tolerated"] = R.frac_basic >= 0.02
R["target"] = R.ref_aa.isin(["D", "E"]) & R.basic_tolerated & ~R.known_site & R.in_safe_region
sel = R[R.target].sort_values(["beta"], ascending=False)
sel.to_csv(D / "mutagenesis_targets.csv", index=False)
pd.set_option("display.width", 220)
print(f"{len(sel)} targets (reference D/E, basic tolerated in >=2 % of bacteria, no documented functional site, YidC restricted to 24-264)\n")
print(sel[["ref", "pos", "ref_aa", "n_aligned", "frac_acidic", "frac_basic", "ph_min_acidic", "ph_min_basic", "beta", "se", "p", "q"]].round(3).to_string(index=False))
print("\nby reference:", sel.groupby("ref").size().to_dict())
print("positions excluded as documented functional sites:", {r: sorted(set(R[(R.ref == r) & R.known_site].pos)) for r in ACC})
print(f"\nper-position associations: {int(msk.sum())} estimated, {int((R.q.dropna() < 0.1).sum())} with q < 0.1, "
      f"{int((R.loc[msk, 'p'] < 0.05).sum())} with p < 0.05 (expected by chance {0.05 * msk.sum():.0f})")
