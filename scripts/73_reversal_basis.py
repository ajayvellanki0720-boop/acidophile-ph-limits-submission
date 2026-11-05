"""Does the experiment's choice of substitution survive the residue decomposition on the panel the headline rests on? -> DATA_DIR/reversal_basis.csv
scripts/36c now runs on both panels. Its hand-picked result was the basis for two design choices that docs/CHARGE_REVERSAL_DESIGN.md and
docs/PREREGISTRATION.md both commit to, and the rule-selected panel changes the support for both.
  1. REVERSAL RATHER THAN NEUTRALISATION. The design replaces D/E with K, not with an uncharged residue. A reversal earns its extra construction cost
     only if the basic term contributes on its own, over and above removing the acidic one. On 7 references K pooled at -0.035 against D+E at +0.041, so
     the basic half was almost as large as the acidic half and reversal looked close to twice as effective as neutralisation. On 32 references K pools at
     -0.010 against D+E at +0.040: the acidic half is unchanged and the basic half is a third of its former size.
  2. D AND E TREATED AS INTERCHANGEABLY ELIGIBLE. scripts/37b selects any position where the E. coli reference carries D or E. On 32 references
     glutamate pools at +0.038 and is positive in 31 of 32, while aspartate pools at +0.010 and is positive in 18 of 32 - a coin flip. The two are not
     the same predictor and the selection criterion treats them as one.
THE QUANTITIES HERE ARE RATIOS, DELIBERATELY. This project has made the same scale error five times (corrections 4, 6, 7, 12 and the per-SD/per-unit
trap in scripts/39a), and the absolute predicted effect is already committed in docs/PREREGISTRATION.md section 3 at the conservative -0.05. Re-deriving
an absolute prediction here would be exactly the post-hoc rescaling that document's section 7 forbids. So nothing below is a new predicted effect: every
figure is a dimensionless ratio between two substitution choices evaluated under the SAME model on the SAME panel, in which the per-position scaling
factor cancels. Flipping one of N exposed periplasmic positions changes the D+E fraction by -1/N and the K fraction by +1/N, so that 1/N divides out of
every ratio and none of them depends on N, on the assay, or on which panel's absolute scaling is believed. The one ratio a reader would most expect - what
the arm's five aspartate positions cost it - is WITHHELD for the same reason, because it can only be formed by mixing coefficients from two separately
fitted univariate models; part 4 records the number it would have given and why that number is wrong.
WHAT THIS DOES NOT DO. It does not say the experiment should be redesigned, and it cannot: whether K is easier to build, better tolerated or more
interpretable than an uncharged substitution is a biochemical question, and a smaller predicted margin for reversal is not an argument for neutralisation
if neutralisation is harder to read. What it does is remove a justification that was resting on 7 references, and record the margin on 32.
Usage: python3 scripts/73_reversal_basis.py   (env DATA_DIR=results/broad2)"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
rows = []

P = {"7-reference hand-picked": pd.read_csv(D / "exposure_decomposition_pooled.csv"),
     "32-reference rule-selected": pd.read_csv(D / "exposure_decomposition_pooled_wide.csv")}

# CORRECTION 15 SUPERSEDES THIS SCRIPT'S RATIOS. Parts 2 and 2b below difference the D+E and lysine coefficients to
# compare a reversal against a neutralisation, and treat both as per percentage point because scripts/36c stores them
# in a column named `beta_per_pt`. They are not: ph.pgls_scan standardises the predictor it tests, so both are per SD
# of their OWN share, and those SDs differ (4.65 against 3.69 points on the rule-selected panel). Differencing them
# mixes two scales, which is the very error the rest of this script was written to avoid - the seventh instance here.
# scripts/74 derives the per-unit coefficients and recomputes every ratio; its figures govern. Kept rather than
# deleted because the ratio printed below is quoted in this project's record as the wrong value.
SUPERSEDED = {"7-reference hand-picked": (2.62, 0.38), "32-reference rule-selected": (1.27, 0.79)}
print("NOTE: parts 2 and 2b below are computed on a mislabelled scale and are superseded by scripts/74 (correction 15).")
print("      Corrected, per unit: the reversal margin is 2.62x on 7 references and 1.27x on 32, and the charge-neutral")
print("      control should recover 38% against 79% - a 40-point gap, not the 25 points these parts report.\n")


def g(df, comp, col="pooled_beta_per_pt"):
    s = df[df.component == comp]
    assert len(s) == 1, f"expected exactly one row for {comp!r}, got {len(s)}"
    return float(s.iloc[0][col])


print("1. the pooled residue coefficients the two panels give (per SD of that component's share - NOT per point,")
print("   despite the column's name; see correction 15 and scripts/74 for the per-point values)")
print(f"   {'component':14s} {'7 refs':>22s} {'32 refs':>22s}")
for comp in ("D", "E", "D+E", "K", "R", "K+R", "net D+E-K-R"):
    cells = []
    for name, df in P.items():
        s = df[df.component == comp].iloc[0]
        cells.append(f"{s.pooled_beta_per_pt:+.3f} ({int(s.n_refs_positive)}/{int(s.n_refs)})")
    print(f"   {comp:14s} {cells[0]:>22s} {cells[1]:>22s}")
    for name, df in P.items():
        s = df[df.component == comp].iloc[0]
        rows.append({"part": "pooled coefficient", "panel": name, "component": comp,
                     "pooled_beta_per_pt": float(s.pooled_beta_per_pt), "pooled_se": float(s.pooled_se),
                     "pooled_p": float(s.pooled_p), "n_refs": int(s.n_refs), "n_refs_positive": int(s.n_refs_positive)})
print("   The sign convention: a POSITIVE coefficient means more of that residue goes with a HIGHER pH_min, i.e. less")
print("   acid tolerance. So acidic residues raise the limit, and removing them or adding basic ones lowers it. Both")
print("   halves of a reversal push the same way, which is the design's rationale for reversing rather than neutralising.")

print("\n2. reversal against neutralisation, as a ratio in which the per-position factor cancels")
print("   Flipping one position from acidic to K changes the D+E fraction by -1/N and the K fraction by +1/N, so the")
print("   predicted shift is (1/N)(b_K - b_DE) for a reversal and (1/N)(-b_DE) for removing the acid alone. N cancels.")
print(f"   {'panel':30s} {'reversal':>10s} {'neutralise':>12s} {'ratio':>8s}")
for name, df in P.items():
    b_de, b_k = g(df, "D+E"), g(df, "K")
    rev, neu = b_k - b_de, -b_de                      # both in units of (pH units) x N, so only their ratio is used
    ratio = rev / neu
    # the identity this rests on, asserted rather than trusted: the reversal effect is the neutralisation effect plus
    # the basic term's own contribution, and nothing here multiplies or divides by a standard deviation
    assert abs(rev - (neu + b_k)) < 1e-12, "reversal must decompose into neutralisation plus the basic term"
    print(f"   {name:30s} {rev:>+10.4f} {neu:>+12.4f} {ratio:>8.2f}x")
    # The reciprocal is the quantity the design's CHARGE-NEUTRAL control tests: how much of the reversal's shift
    # removing the acid alone should recover. It is stored rather than left to be computed in prose from the ratio,
    # because a derived quantity quoted anywhere has to be a row or it is unguarded (correction 12).
    rows.append({"part": "reversal vs neutralisation SUPERSEDED", "panel": name, "beta_DE": b_de, "beta_K": b_k,
                 "effect_reversal_times_N": float(rev), "effect_neutralise_times_N": float(neu),
                 "ratio_reversal_over_neutralise": float(ratio),
                 "share_neutral_control_of_reversal": float(neu / rev),
                 "ratio_corrected_scripts74": SUPERSEDED[name][0],
                 "share_corrected_scripts74": SUPERSEDED[name][1],
                 "note": "computed on per-SD coefficients mislabelled per point; correction 15, see scripts/74"})
r7 = (g(P["7-reference hand-picked"], "K") - g(P["7-reference hand-picked"], "D+E")) / -g(P["7-reference hand-picked"], "D+E")
r32 = (g(P["32-reference rule-selected"], "K") - g(P["32-reference rule-selected"], "D+E")) / -g(P["32-reference rule-selected"], "D+E")
print(f"   So the margin reversal buys over neutralisation falls from {r7:.2f}x to {r32:.2f}x. The acidic half is intact")
print(f"   ({g(P['7-reference hand-picked'], 'D+E'):+.3f} against {g(P['32-reference rule-selected'], 'D+E'):+.3f}); it is the basic half that does not replicate.")
rows.append({"part": "verdict reversal", "ratio_7ref": float(r7), "ratio_32ref": float(r32),
             "note": "margin of reversal over neutralisation, scale-free; falls because the basic term shrinks"})
print("\n2b. this is a PREDICTION the experiment can already test, with a control it is going to build anyway")
print("   docs/CHARGE_REVERSAL_DESIGN.md section 6 includes a charge-neutral control (D->N, E->Q at the same positions)")
print("   whose stated expectation is only that it lands 'intermediate'. The two panels make that vague expectation")
print("   sharp, and they disagree: the share of the reversal's shift that neutralisation alone should recover is")
for name, df in P.items():
    b_de, b_k = g(df, "D+E"), g(df, "K")
    share = -b_de / (b_k - b_de)
    print(f"      {name:30s} {share:.0%}")
    rows.append({"part": "neutral control prediction", "panel": name, "beta_DE": b_de, "beta_K": b_k,
                 "share_neutral_control_of_reversal": float(share)})
s7 = -g(P["7-reference hand-picked"], "D+E") / (g(P["7-reference hand-picked"], "K") - g(P["7-reference hand-picked"], "D+E"))
s32 = -g(P["32-reference rule-selected"], "D+E") / (g(P["32-reference rule-selected"], "K") - g(P["32-reference rule-selected"], "D+E"))
print(f"   {s7:.0%} against {s32:.0%} is a {abs(s32 - s7):.0%}-point gap in a quantity the control measures directly, so the")
print("   control stops being a mere sanity check and becomes the arbiter of which panel's residue decomposition is")
print("   right. Stated here, before any strain exists, so it cannot be fitted to whatever the control turns out to do.")
print(f"   The honest caveat: distinguishing {s7:.0%} from {s32:.0%} of a shift that is itself predicted at -0.05 pH units needs")
print("   resolution the confirmation assay does not have, so this is a prediction the POOLED assay might test if its")
print("   Gate 2 slope is ever measured, and not one the primary endpoint can settle.")
rows.append({"part": "verdict neutral control", "share_7ref": float(s7), "share_32ref": float(s32),
             "gap_percentage_points": float(abs(s32 - s7) * 100),
             "note": "prediction recorded before construction; needs pooled-assay resolution, not the primary endpoint"})

print("\n3. the same comparison on the net composite, which is a DIFFERENT model and must not be mixed with part 2")
print("   The headline predictor is the single net term D+E-K-R. Under that one model a reversal moves net by -2/N and")
print("   a neutralisation by -1/N, so the ratio is exactly 2 by construction and carries no information about K at all.")
for name, df in P.items():
    b_net = g(df, "net D+E-K-R")
    rows.append({"part": "net basis", "panel": name, "beta_net": b_net,
                 "ratio_reversal_over_neutralise": 2.0,
                 "note": "ratio is 2 by construction under the single net predictor; not comparable with the separate-term ratio"})
    print(f"   {name:30s} net {b_net:+.3f}; ratio fixed at 2.00x by the predictor's definition")
print("   This is why part 2 uses the separate terms: only a model that gives D+E and K their own coefficients can say")
print("   whether the basic half contributes. Quoting the net model's 2x as evidence for reversal would be circular.")

print("\n4. the primary arm's own composition, against a criterion that treats D and E alike")
T = pd.read_csv(D / "mutagenesis_targets.csv")
T = T[T.target]
arm = T[T.ref.isin(["NuoL", "NuoN"])]
assert len(arm) == 14, f"the primary arm is 14 positions in docs/PREREGISTRATION.md section 2; found {len(arm)}"
nE, nD = int((arm.ref_aa == "E").sum()), int((arm.ref_aa == "D").sum())
print(f"   NuoL + NuoN combined variant: {len(arm)} positions, {nE} glutamate and {nD} aspartate")
for name, df in P.items():
    b_e, b_d, b_de = g(df, "E"), g(df, "D"), g(df, "D+E")
    rows.append({"part": "arm composition", "panel": name, "n_positions": int(len(arm)), "n_glutamate": nE,
                 "n_aspartate": nD, "beta_E": b_e, "beta_D": b_d, "beta_DE": b_de})
# NO RATIO IS REPORTED HERE, and the reason is the whole point of the project's correction history. The obvious
# quantity - weighting each position by its own residue, nE*b_E + nD*b_D, against treating all 14 as D+E, 14*b_DE -
# mixes coefficients from DIFFERENT models. scripts/36c fits each component on its own, so b_E comes from a model in
# which aspartate variation is unmodelled and b_DE from one in which the two are summed into a single predictor;
# b_E and b_D are not additive contributions that sum to b_DE, and nothing licenses adding them. Computing it anyway
# gave 0.63x on the hand-picked panel and 0.71x on the rule-selected one - below 1 on BOTH, which is the tell, since
# a genuine consequence of aspartate's failure to replicate would show up on the wide panel and not on the narrow one.
# This is the sixth instance of this project's recurring error and the first caught before the figure was published.
# The legitimate version needs D and E entered SIMULTANEOUSLY, which scripts/36c does not currently fit; until it
# does, the finding below is stated qualitatively, which is all the available fits support.
rows.append({"part": "ratio withheld", "panel": "both",
             "note": "per-position residue weighting not reported: b_D and b_E come from separate univariate fits and "
                     "are not additive against b_DE from the summed-predictor model; needs a simultaneous D,E fit"})
e_share = nE / len(arm)
print(f"   The arm is {e_share:.0%} glutamate, which is the favourable direction: on 32 references E is the component that")
print(f"   replicates (31/32) and D is the one that does not (18/32). That is an accident of selection rather than a")
print(f"   design choice - scripts/37b ranks positions on each one's own association and never on which residue it is -")
print(f"   so the criterion's indifference between D and E remains a real weakness that this particular arm escaped.")
print(f"   HOW MUCH it costs is not quantified here; see the comment above on why that ratio is withheld. The {nD}")
print(f"   aspartate positions are listed so a reviewer can see which part of the arm rests on the weaker component.")
print(f"   {'ref':>6} {'pos':>5} {'aa':>3}   (the aspartate positions, whose residue pools at 18/32)")
for x in arm[arm.ref_aa == "D"].itertuples():
    print(f"   {x.ref:>6} {int(x.pos):>5} {x.ref_aa:>3}")
    rows.append({"part": "aspartate position", "panel": "both", "ref": str(x.ref), "pos": int(x.pos), "ref_aa": "D"})

R = pd.DataFrame(rows)
R.to_csv(D / "reversal_basis.csv", index=False)
print("\nwrote", D / "reversal_basis.csv")
print("\nReading. The acidic half of the design's rationale is confirmed and sharpened by the wider panel: D+E pools at")
print("essentially the same value with a p five orders of magnitude smaller, and the carrier is specifically glutamate.")
print("The basic half does not replicate - K falls from -0.035 to -0.010 and from 5/6 to 11/32 - so the margin that")
print("reversal buys over simply removing the acid narrows sharply (2.62x to 1.27x on the corrected per-unit scale of")
print("scripts/74; the 1.86x to 1.26x printed above is on the mislabelled scale). That is not a reason")
print("to change the substitution, which is a construction decision this analysis cannot make, but the margin should no")
print("longer be quoted from the hand-picked panel, and the eligibility criterion's indifference between D and E is a")
print("real weakness that this arm happens to have escaped, being 64% glutamate by accident of ranking.")
print("One quantity a reader might expect is deliberately absent: how much the five aspartate positions cost the arm.")
print("Producing it requires a model with D and E entered simultaneously, and the available fits are univariate per")
print("component, so the arithmetic would mix models. The comment in part 4 records the figure it would have given and")
print("why it is wrong, because a withheld number is easier to re-derive by accident than a refuted one.")
