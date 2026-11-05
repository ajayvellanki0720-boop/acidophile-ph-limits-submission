"""Is the exposed-periplasmic signal distributed over many positions on the panel that carries the headline? -> DATA_DIR/distribution_wide.csv and DATA_DIR/distribution_wide_by_reference.csv
scripts/36c's position-distribution test was the last claim in docs/SYNTHESIS.md section 3 still measured only on the 7-reference hand-picked panel, where it read as
settled: both random halves of a reference's exposed periplasmic positions came out positive in 92% of splits for YidC, 95% for NuoN and 85% for NuoL, against 30% for
cytochrome bo3. Run on the 32 rule-selected references the headline number is much worse - the median reference agrees across halves in only 45% of splits and just 13 of 32
reach 80% - and the question is whether that refutes the claim or reflects something else.
IT REFLECTS SOMETHING ELSE, AND THE TEST THAT SHOWS SO IS NOT CIRCULAR. The statistic is conditional on there being a signal to distribute, and most of the 32 references
have none of their own: scripts/54 established that a per-reference estimate is largely a measure of how cleanly that reference can be scored, and only 27 of 32 are even
positive. Among the 10 references individually significant at p<0.05 the agreement is 100% median, 10 of 10 at or above 80%; among the other 22 it is 35% median, 3 of 22.
WHY THE CORRELATION WITH EFFECT SIZE IS NOT MECHANICAL, which is the objection to answer. A large coefficient does not force cross-half agreement. If one position carried a
reference's whole effect, each random half would contain it about half the time, so the two halves would agree in sign in roughly 50% of splits NO MATTER how large the
coefficient was. Agreement near 100% therefore means the effect survives being cut in half, which is what "distributed" asserts; agreement near 50% with a large coefficient
is the concentrated signature; agreement near 25-40% is no signal either way. So the statistic discriminates between concentration and distribution independently of size,
and the correct reading of a reference at 35% with a near-zero coefficient is "nothing to distribute", not "concentrated".
THE STATISTIC IS SIGN-NEUTRAL HERE, unlike 36c's. 36c reports the share of splits in which both halves are POSITIVE, which sends the five negative-coefficient references to
near zero because their halves agree negatively - that measures sign, not distribution. This uses the share in which the halves AGREE, either way.
WHAT THIS COSTS THE CLAIM, stated plainly. "The signal is distributed rather than carried by a few sites" has been asserted here without qualification. It should read: where
a reference has a signal, that signal is distributed across its positions; across the panel as a whole most references have no individual signal to speak of. One reference,
P62517, shows the concentrated signature - a coefficient of reasonable size with halves agreeing only about half the time - and is named rather than averaged away.
Usage: python3 scripts/78_distribution_wide.py   (env DATA_DIR=results/broad2)"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
BIG, MID = 0.03, 0.65          # the concentrated signature: a coefficient this large with agreement this low
rows = []

P = pd.read_csv(D / "exposure_position_splits_wide.csv")
E = pd.read_csv(D / "exposure_decomposition_wide.csv")
net = (E[E.component == "net D+E-K-R"][["ref", "beta_per_pt", "sd_pts", "p"]]
       .rename(columns={"beta_per_pt": "beta_ref", "p": "p_ref"}))
g = P.dropna(subset=["beta_a", "beta_b"]).copy()
g["same"] = np.sign(g.beta_a) == np.sign(g.beta_b)
S = g.groupby("ref")["same"].mean().rename("share_same_sign").reset_index()
S = S.merge(g.groupby("ref")["positions"].first().reset_index(), on="ref")
M = S.merge(net, on="ref")
assert len(M) >= 30, f"expected the 32-reference panel; matched {len(M)}"

print(f"1. cross-half agreement across {len(M)} references (share of 40 random splits in which both halves agree in sign)")
print(f"   median {M.share_same_sign.median():.0%}, range {M.share_same_sign.min():.0%} to {M.share_same_sign.max():.0%}, "
      f"{int((M.share_same_sign >= 0.8).sum())} of {len(M)} at 80% or above")
print("   On the hand-picked panel of 7 the comparable figures were 92%, 95% and 85% for three of four references,")
print("   which is why the claim has been stated here without qualification. It cannot be, on this panel.")
rows.append({"part": "overall", "n_refs": int(len(M)), "median_share": float(M.share_same_sign.median()),
             "min_share": float(M.share_same_sign.min()), "max_share": float(M.share_same_sign.max()),
             "n_at_80pct": int((M.share_same_sign >= 0.8).sum())})

print("\n2. conditional on the reference having a signal of its own, which is what the statistic requires")
print(f"   {'subset':34s} {'n':>3s} {'median agreement':>17s} {'at >=80%':>9s}")
for lab, sub in (("individually significant, p<0.05", M[M.p_ref < 0.05]),
                 ("not individually significant", M[M.p_ref >= 0.05]),
                 ("coefficient positive", M[M.beta_ref > 0]),
                 ("coefficient zero or negative", M[M.beta_ref <= 0])):
    if not len(sub):
        continue
    print(f"   {lab:34s} {len(sub):>3d} {sub.share_same_sign.median():>17.0%} "
          f"{str(int((sub.share_same_sign >= 0.8).sum())) + '/' + str(len(sub)):>9s}")
    rows.append({"part": "by subset", "subset": lab, "n_refs": int(len(sub)),
                 "median_share": float(sub.share_same_sign.median()),
                 "n_at_80pct": int((sub.share_same_sign >= 0.8).sum())})
r_abs = float(M.share_same_sign.corr(M.beta_ref.abs()))
r_pos = float(M.share_same_sign.corr(M.positions))
print(f"   agreement against |coefficient|: r {r_abs:+.3f}; against class size: r {r_pos:+.3f}")
print("   So it tracks whether the reference has an effect, not how many positions it has.")
rows.append({"part": "correlations", "r_agreement_vs_abs_beta": r_abs, "r_agreement_vs_positions": r_pos})

print("\n3. the objection: is that correlation mechanical? No, and here is the reason it is not")
print("   If ONE position carried a reference's whole effect, each random half would contain it about half the")
print("   time, so the halves would agree in roughly 50% of splits however large the coefficient. Agreement near")
print("   100% therefore means the effect survives being cut in half, which is what 'distributed' asserts.")
print("   Agreement near 50% WITH a sizeable coefficient is the concentrated signature, and it is distinguishable:")
conc = M[(M.beta_ref.abs() > BIG) & (M.share_same_sign < MID)]
print(f"   references with |coefficient| > {BIG} and agreement below {MID:.0%}: {len(conc)}")
for x in conc.itertuples():
    print(f"      {x.ref:>8} beta {x.beta_ref:+.4f}, agreement {x.share_same_sign:.0%}, {int(x.positions)} positions")
    rows.append({"part": "concentrated signature", "ref": str(x.ref), "beta_ref": float(x.beta_ref),
                 "share_same_sign": float(x.share_same_sign), "positions": int(x.positions)})
if not len(conc):
    print("      (none)")
print("   A reference in that band is the case the test exists to catch, and it is reported rather than averaged")
print("   into the median.")

M.sort_values("share_same_sign", ascending=False).to_csv(D / "distribution_wide_by_reference.csv", index=False)
pd.DataFrame(rows).to_csv(D / "distribution_wide.csv", index=False)
print("\nwrote", D / "distribution_wide.csv", "and", D / "distribution_wide_by_reference.csv")
print("\nReading. The unqualified claim has to go. 'The signal is distributed rather than carried by a few sites' is")
print("supported where there is a signal - 10 of 10 individually significant references agree across halves at 80% or")
print("better, median 100% - and says nothing about the 22 references with no individual effect, which agree at 35%")
print("median because there is nothing to distribute. The panel-wide median of 45% is therefore not evidence against")
print("distribution; it is a restatement of scripts/54's finding that most per-reference estimates are scoring quality")
print("rather than biology. The honest form of the claim is conditional, and one reference looks concentrated.")
