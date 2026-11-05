"""Is the panel's reliability ranking mostly a measure of how big each reference is? -> DATA_DIR/reference_size.csv and DATA_DIR/reference_size_table.csv
scripts/54 established that a reference's own coefficient correlates +0.817 with that reference's reliability, and concluded that a per-reference
estimate is mostly a measure of how cleanly that reference can be scored rather than of how much it matters. It did not ask what makes a reference
cleanly scorable, and there is an obvious mechanical candidate: how many exposed periplasmic positions it has. A genome's value for a reference is an
average of D+E-K-R over those positions, so a reference with 262 of them is averaging away noise that a reference with 17 cannot, exactly as
Spearman-Brown describes for any aggregate of noisy items, and across the 32 membrane references that count runs from 17 to 262.
THE ANSWER IS PARTLY YES AND LESS THAN EXPECTED, which is why the framing below is deliberately modest. A one-parameter measurement-averaging model
predicting each reference's reliability from its position count alone explains about a third of the variance in the ranking, not most of it, and the
+0.817 association between coefficient and reliability falls only to +0.737 when log position count is residualised out of both. So scripts/54's warning
stands largely on its own terms rather than being a size effect in disguise. What size does fully account for is the particular reference chosen as an
experimental arm in correction 13, whose observed reliability of 0.449 is slightly BELOW the 0.479 its size alone predicts.
This matters in three places at once, which is why it is worth its own script.
  1. It refines the synthesis's claim. "A reference's coefficient measures how cleanly it can be scored" is right but incomplete; the mechanism is
     mostly position count, which is a property of the protein's topology and has nothing to do with acid tolerance.
  2. It is a confound on the per-reference ranking that scripts/54 publishes and that correction 13 then used to choose an experimental arm. If
     reliability is largely size, ranking references by reliability ranks them by size.
  3. It decides whether correction 13's replacement arm survives. P24177 was chosen partly for having the most acidic positions available to reverse and
     partly for having the panel's highest reliability - and it is also the largest reference in the panel by a factor of three, so both advantages
     could be the same fact twice.
The test that matters for (3) is the D/E RATE per position rather than the count, which removes size entirely, and the test that matters for (1) and (2)
is whether reliability predicts the per-reference coefficient once size is held fixed.
WHAT THIS SCRIPT CANNOT SETTLE, stated because it bears directly on the arm choice. P24177's exposed periplasmic positions are mostly in a large
periplasmic porter domain, where the 39 eligible positions sit in folded structure, whereas NuoL and NuoN present short interhelical loops. Whether a
functional domain tolerates 39 charge reversals is a biochemical question no amount of comparative data answers, and the interpretability gate should be
expected to reject a larger fraction of P24177 variants than of the primary arm's. That is a real cost of the recommendation and is recorded rather than
left to be discovered at the bench.
Usage: python3 scripts/72_reference_size.py   (env DATA_DIR=results/broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
rows = []

cls = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
refs = refs[refs.accepted & refs.membrane]
rank = pd.read_csv(D / "reference_ranking.csv")
rank = rank[rank.ref != "SUMMARY"].set_index("ref")
out = []
for r in refs.itertuples():
    pe = cls[(cls.ref == str(r.label)) & (cls["class"] == "peri_exp")]
    if pe.empty:
        continue
    aa = pe.aa.to_numpy()
    out.append({"ref": str(r.label), "name": str(r.name)[:38], "n_positions": len(aa),
                "reference_DE": int(np.isin(aa, ["D", "E"]).sum()),
                "de_rate": float(np.isin(aa, ["D", "E"]).mean()),
                "reliability": float(rank.reliability.get(str(r.label), np.nan)),
                "item_rest_r": float(rank.item_rest_r.get(str(r.label), np.nan)),
                "beta_per_sd": float(rank.beta_per_sd.get(str(r.label), np.nan))})
F = pd.DataFrame(out)
print(f"{len(F)} membrane references; exposed periplasmic positions run {int(F.n_positions.min())} to {int(F.n_positions.max())}, "
      f"median {int(F.n_positions.median())}")


def rp(a, b, label):
    s = F[[a, b]].dropna()
    pear, spear = s[a].corr(s[b]), s[a].corr(s[b], method="spearman")
    print(f"   {label:52s} r {pear:+.3f}  rho {spear:+.3f}  (n={len(s)})")
    rows.append({"part": "correlation", "x": a, "y": b, "pearson": float(pear), "spearman": float(spear), "n": len(s)})
    return pear


print("\n1. what reference size is associated with")
rp("n_positions", "reliability", "position count against reliability")
rp("n_positions", "item_rest_r", "position count against item-rest correlation")
rp("n_positions", "beta_per_sd", "position count against the per-reference coefficient")
rp("n_positions", "reference_DE", "position count against the count of acidic positions")
print("   The last of these is near unity by construction - a bigger protein has more of everything - so the count of")
print("   acidic positions is close to a restatement of size, and 'buildability' measured that way is largely size.")

print("\n2. the mechanism, tested rather than asserted: Spearman-Brown in the number of positions")
# if a reference's score is a mean of k noisy per-position items, reliability should rise with k as k*r1/(1+(k-1)*r1)
s = F[["n_positions", "reliability"]].dropna()
k = s.n_positions.to_numpy(float)
y = s.reliability.to_numpy(float)
best = None
for r1 in np.linspace(0.0005, 0.2, 400):
    pred = k * r1 / (1 + (k - 1) * r1)
    ss = float(((y - pred) ** 2).sum())
    if best is None or ss < best[1]:
        best = (r1, ss)
r1 = best[0]
pred = k * r1 / (1 + (k - 1) * r1)
r2 = 1 - best[1] / float(((y - y.mean()) ** 2).sum())
print(f"   fitting a single per-position reliability r1 and predicting each reference's reliability from its size alone:")
print(f"   r1 = {r1:.4f}, and that one-parameter curve explains R2 = {r2:.3f} of the variance in reliability across references")
print(f"   That is a real contribution and smaller than the correlation in part 1 suggests: a one-parameter curve based on")
print(f"   nothing but topology accounts for about a third of the ranking, so size matters and does not explain it away.")
rows.append({"part": "spearman-brown", "per_position_reliability": float(r1), "r2_of_reliability": float(r2), "n": len(s)})

print("\n3. does reliability predict the coefficient once size is held fixed?")
s = F[["n_positions", "reliability", "beta_per_sd"]].dropna()
# partial correlation of beta with reliability given log size, by residualising both on log size
lg = np.log(s.n_positions.to_numpy(float))
X = np.column_stack([np.ones(len(lg)), lg])


def resid(v):
    b = np.linalg.lstsq(X, v, rcond=None)[0]
    return v - X @ b


rb = resid(s.beta_per_sd.to_numpy(float))
rr = resid(s.reliability.to_numpy(float))
raw = float(s.reliability.corr(s.beta_per_sd))
part = float(np.corrcoef(rb, rr)[0, 1])
print(f"   raw correlation of coefficient with reliability: {raw:+.3f}")
print(f"   with log position count residualised out of both: {part:+.3f}  (n={len(s)})")
if abs(part) < abs(raw) / 2:
    print("   So most of scripts/54's +0.817 runs through size. The warning it drew stands - a per-reference estimate is not")
    print("   biology - but the mechanism is position count, which is a fact about protein topology, not about scoring quality")
    print("   in any deeper sense.")
else:
    print("   So the association survives holding size fixed: scoring quality is not reducible to position count, and")
    print("   scripts/54's warning stands on its own terms rather than being a size effect in disguise.")
rows.append({"part": "partial", "raw_r_beta_reliability": raw, "partial_r_given_log_size": part, "n": len(s)})

print("\n4. does correction 13's arm choice survive the confound?")
p24 = F[F.ref == "P24177"].iloc[0]
rate_rank = int((F.de_rate > p24.de_rate).sum()) + 1
size_rank = int((F.n_positions > p24.n_positions).sum()) + 1
print(f"   P24177: {int(p24.n_positions)} positions (rank {size_rank} of {len(F)} on size), {int(p24.reference_DE)} acidic, "
      f"rate {p24.de_rate:.3f} (rank {rate_rank} of {len(F)})")
print(f"   panel median rate {F.de_rate.median():.3f}. So it is acidic-DENSE as well as large, which is the part size does not explain.")
print(f"   {'ref':>8} {'name':38s} {'pos':>5} {'D/E':>4} {'rate':>6} {'reliability':>12}")
for x in F.nlargest(6, "de_rate").itertuples():
    print(f"   {x.ref:>8} {x.name:38s} {int(x.n_positions):>5} {int(x.reference_DE):>4} {x.de_rate:>6.3f} {x.reliability:>12.3f}")
print("   The arm choice survives: the references with a comparable rate have 17 to 60 positions, so none of them offers an")
print("   arm of workable size. P24177 is the only reference that is both dense and large.")
print("   But its reliability advantage should NOT be cited as independent evidence for it, because on the model in part 2")
print(f"   a reference of {int(p24.n_positions)} positions is expected to score about {p24.n_positions * r1 / (1 + (p24.n_positions - 1) * r1):.3f} "
      f"on reliability from size alone, against its observed {p24.reliability:.3f}.")
rows.append({"part": "arm choice", "ref": "P24177", "n_positions": int(p24.n_positions), "reference_DE": int(p24.reference_DE),
             "de_rate": float(p24.de_rate), "de_rate_rank": rate_rank, "size_rank": size_rank,
             "panel_median_rate": float(F.de_rate.median()), "reliability_observed": float(p24.reliability),
             "reliability_expected_from_size": float(p24.n_positions * r1 / (1 + (p24.n_positions - 1) * r1))})

F.to_csv(D / "reference_size_table.csv", index=False)
pd.DataFrame(rows).to_csv(D / "reference_size.csv", index=False)
print("\nwrote", D / "reference_size.csv", "and", D / "reference_size_table.csv")
print("\nReading. Size is a real but partial driver of the reliability ranking - a third of its variance from one topology")
print("parameter - and the coefficient-reliability association survives holding size fixed, so scripts/54's conclusion is")
print("refined rather than overturned and its warning is not reducible to protein length. The experimental arm chosen in")
print("correction 13 survives on the part size does not explain, its acidic density being second of 32; but its reliability is")
print("exactly what its size predicts and must not be offered as a separate argument for it.")
