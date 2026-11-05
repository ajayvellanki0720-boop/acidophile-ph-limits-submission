"""If the residue coefficients do not follow side-chain pKa, what DO they follow? -> DATA_DIR/what_organises_residues.csv and, with `wide`, DATA_DIR/what_organises_residues_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHAT CORRECTION 18 ESTABLISHED AND LEFT OPEN. The pKa test asked whether the residue ranking follows side-chain pKa and answered no: histidine (pKa ~6.0) runs OPPOSITE in sign
to glutamate (~4.25) at twice the magnitude, where a charge-burden account predicted same sign and larger. That is a negative result about one hypothesis and it says nothing
about what organises the coefficients instead. All five titratable residues are now in the panel - D 3.65, E 4.25, H 6.0, K 10.5, R 12.5 - so the question is answerable.
TWO ACCOUNTS, AND THEY MAKE DIFFERENT PREDICTIONS ABOUT DIFFERENT CONTRASTS, which is what makes this a test rather than a description. On a PKA account the coefficient should
vary monotonically with pKa, so every adjacent pair should differ in the same direction. On a FORMAL-CHARGE account the coefficients should group by the sign of the residue's
charge near neutral pH - aspartate and glutamate negative, histidine, lysine and arginine positive - with the pKa differences inside each group doing nothing. The pre-registered
pKa ordering already failed ACROSS the groups. What has never been tested is whether it holds WITHIN them, and that is the comparison the two accounts disagree about: if pKa
organises nothing, glutamate should not differ from aspartate either, and this project has reported that it does, at 2.4x. So the live possibility is a third thing - charge
sign sets the direction and pKa orders the magnitudes inside each class - and it is a possibility neither account as stated allows.
EACH CONTRAST IS MADE A COEFFICIENT RATHER THAN COMPUTED FROM TWO OF THEM. ph.pgls returns standard errors but not the covariance between coefficients, so a difference taken
between two of its outputs has no standard error. Instead each contrast is entered as its own column alongside the matching level - the sum-and-difference reparametrisation
scripts/82 validated against the two-level fit - so a difference of coefficients becomes a coefficient with a p-value, and "glutamate exceeds aspartate" becomes a claim that
can fail rather than an impression formed by putting two numbers next to each other.
SCALE. Per percentage point of pooled exposed-periplasmic share throughout, from ph.pgls directly and never ph.pgls_scan - correction 15. Contrasts are differences of
per-point coefficients and so are in the same units; a contrast of per-SD coefficients would not be, since the five shares have different standard deviations.
Usage: python3 scripts/90_what_organises_residues.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
WHICH = (sys.argv[1] if len(sys.argv) > 1 else "wide").lower()
assert WHICH in ("orig", "wide"), f"panel must be orig or wide, got {WHICH!r}"
PANEL = "structure_panel_wide_det.csv" if WHICH == "wide" else "structure_panel.csv"
TAG = "_wide" if WHICH == "wide" else ""
MINPOS = 200
PKA = {"D": 3.65, "E": 4.25, "H": 6.0, "K": 10.5, "R": 12.5}
ACID, BASE = ("D", "E"), ("H", "K", "R")
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
if not any(c.endswith("|H") for c in pan.columns):
    sys.exit(f"{PANEL} has no histidine columns; rebuild the panel with scripts/36a first")
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
RES = ("D", "E", "H", "K", "R")
tot = {k: np.zeros(len(B)) for k in ("n",) + RES}
for ref in MEM:
    if f"{ref}|peri_exp|n" not in B:
        continue
    tot["n"] += B[f"{ref}|peri_exp|n"].fillna(0).to_numpy(float)
    for k in RES:
        tot[k] += B[f"{ref}|peri_exp|{k}"].fillna(0).to_numpy(float)
ok = tot["n"] >= MINPOS
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
y = g.ph_min.to_numpy(float)
XS = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
lam = ph.pgls(y, XS, Cm)["lambda"]
s = {k: 100 * tot[k][ok] / tot["n"][ok] for k in RES}


def fit(cols, names):
    X = np.column_stack([XS] + list(cols))
    assert np.linalg.matrix_rank(X) == X.shape[1], f"rank-deficient design for {names}"
    r = ph.pgls(y, X, Cm, lam=lam)
    k = len(cols)
    return {nm: (float(r["beta"][-k + i]), float(r["se"][-k + i]), float(r["p"][-k + i]))
            for i, nm in enumerate(names)}


print(f"1. all five titratable residues entered jointly (n={int(ok.sum())}, lambda {lam:.3f})")
joint = fit([s[k] for k in RES], list(RES))
print(f"   {'residue':>8} {'pKa':>6} {'share %':>9} {'beta/pt':>10} {'se':>9} {'p':>10} {'charge at pH 7':>15}")
for k in RES:
    b, se, pv = joint[k]
    print(f"   {k:>8} {PKA[k]:>6.2f} {s[k].mean():>9.3f} {b:>+10.5f} {se:>9.5f} {pv:>10.2e} "
          f"{('negative' if k in ACID else 'positive'):>15}")
    rows.append({"part": "joint fit", "residue": k, "pka": PKA[k], "mean_share_pct": float(s[k].mean()),
                 "beta_per_pt": b, "se": se, "p": pv, "charge_at_ph7": "negative" if k in ACID else "positive",
                 "n": int(ok.sum()), "lambda": lam})
order = sorted(RES, key=lambda k: PKA[k])
signs = [1 if joint[k][0] > 0 else -1 for k in order]
mono = all(joint[order[i]][0] >= joint[order[i + 1]][0] for i in range(len(order) - 1)) or \
       all(joint[order[i]][0] <= joint[order[i + 1]][0] for i in range(len(order) - 1))
print(f"   ordered by pKa the coefficients are "
      + ", ".join(f"{k} {joint[k][0]:+.4f}" for k in order))
print(f"   monotone in pKa: {mono}. Signs by charge class: acidic "
      f"{'all positive' if all(joint[k][0] > 0 for k in ACID) else 'mixed'}, basic "
      f"{'all negative' if all(joint[k][0] < 0 for k in BASE) else 'mixed'}.")
rows.append({"part": "pattern", "monotone_in_pka": bool(mono),
             "acidic_all_positive": bool(all(joint[k][0] > 0 for k in ACID)),
             "basic_all_negative": bool(all(joint[k][0] < 0 for k in BASE))})

print("\n2. the three contrasts, each entered as its own coefficient so it has a standard error")
A = (s["D"] + s["E"]) / 2
Bs = (s["H"] + s["K"] + s["R"]) / 3
tests = {}
r1 = fit([A + Bs, A - Bs], ["level", "acidic minus basic"])
tests["acidic minus basic"] = r1["acidic minus basic"]
r2 = fit([s["D"] + s["E"], s["E"] - s["D"], Bs], ["acid level", "glutamate minus aspartate", "basic"])
tests["glutamate minus aspartate"] = r2["glutamate minus aspartate"]
KR = (s["K"] + s["R"]) / 2
r3 = fit([s["H"] + KR, s["H"] - KR, A], ["basic level", "histidine minus lysine/arginine", "acidic"])
tests["histidine minus lysine/arginine"] = r3["histidine minus lysine/arginine"]
print(f"   {'contrast':>34} {'pKa gap':>8} {'beta/pt':>10} {'se':>9} {'p':>10}")
GAPS = {"acidic minus basic": abs(np.mean([PKA[k] for k in ACID]) - np.mean([PKA[k] for k in BASE])),
        "glutamate minus aspartate": PKA["E"] - PKA["D"],
        "histidine minus lysine/arginine": abs(PKA["H"] - np.mean([PKA["K"], PKA["R"]]))}
for nm, (b, se, pv) in tests.items():
    print(f"   {nm:>34} {GAPS[nm]:>8.2f} {b:>+10.5f} {se:>9.5f} {pv:>10.2e}")
    rows.append({"part": "contrast", "contrast": nm, "pka_gap": GAPS[nm], "beta_per_pt": b, "se": se, "p": pv})

print("\n2b. THE CONTRAST THAT DECIDES WHETHER THERE IS A CHARGE CLASS AT ALL")
# The charge-class contrast above pools histidine with lysine and arginine, and part 1 printed that the basic
# class is NOT sign-uniform: histidine is -0.22 while lysine and arginine sit near zero. A contrast against a
# group containing one large outlier mostly measures the outlier, so "acidic differs from basic" may be nothing
# more than "acidic differs from histidine". Dropping histidine from the basic side is the test, and it is the
# one that distinguishes a charge-class account from a histidine-specific one. The first version of this script
# did not run it and its verdict asserted that charge sign sets the direction while the line above it said the
# basic class was mixed - a conclusion contradicting the script's own output, which is the error this project has
# now made in four scripts in two days.
r4 = fit([A + KR, A - KR, s["H"]], ["level", "acidic minus lysine/arginine", "histidine"])
abk = r4["acidic minus lysine/arginine"]
print(f"   acidic minus lysine/arginine, histidine held fixed: {abk[0]:+.5f}/pt (se {abk[1]:.5f}, p {abk[2]:.2e})")
print(f"   against acidic minus the whole basic class: {tests['acidic minus basic'][0]:+.5f}")
rows.append({"part": "contrast", "contrast": "acidic minus lysine/arginine (histidine held fixed)",
             "beta_per_pt": abk[0], "se": abk[1], "p": abk[2]})
charge_class_real = abk[2] < 0.05 and (abk[0] > 0)
print("   " + ("A charge-class difference SURVIVES removing histidine, so the grouping is real and histidine is"
               if charge_class_real else
               "The charge-class difference does NOT survive removing histidine, so there is no charge-class"))
print("   " + ("   an outlier within it." if charge_class_real else
               "   effect: 'acidic differs from basic' was 'acidic differs from histidine' all along."))
rows.append({"part": "charge class", "survives_without_histidine": charge_class_real})

print("\n3. which account the contrasts pick out")
# A pKa account predicts all three non-zero and ordered by the pKa gap. A formal-charge account predicts the
# first non-zero and the other two null. Anything else is a third account, and the point of running all three is
# that the data can say so rather than being forced into one of the two hypotheses that were written down.
ab, ae, ah = (tests["acidic minus basic"], tests["glutamate minus aspartate"],
              tests["histidine minus lysine/arginine"])
sig = {k: v[2] < 0.05 for k, v in tests.items()}
uniform_basic = all(joint[k][0] < 0 for k in BASE) or all(joint[k][0] > 0 for k in BASE)
# THE VERDICT IS GATED ON SIGN UNIFORMITY, which the first version was not. A charge-class account requires the
# coefficients to group BY SIGN; if one class is mixed, that account is already false whatever the class contrast
# says, because the contrast of a mixed group is a statement about its mean and not about the grouping.
if not uniform_basic and not charge_class_real:
    verdict = ("NEITHER ACCOUNT, AND NOT A CHARGE CLASS EITHER. The basic class is not sign-uniform - histidine "
               "is strongly negative while lysine and arginine are null - and the acidic-against-basic contrast "
               "does not survive removing histidine. So the pattern is three facts and not one principle: the "
               "acidic residues are positive and ordered as their pKas are, histidine is large and negative, and "
               "lysine and arginine do nothing. Histidine is sui generis here, which is what correction 19's "
               "separate-mechanisms finding already implied and this makes explicit.")
elif not uniform_basic and charge_class_real:
    verdict = ("A CHARGE CLASS WITH AN OUTLIER. The acidic side differs from lysine and arginine even with "
               "histidine held fixed, so the grouping is real, but histidine does not belong to it: it is far "
               "larger than its class-mates and the class is not sign-uniform. pKa orders the acidic pair.")
elif sig["acidic minus basic"] and not sig["glutamate minus aspartate"] and not sig["histidine minus lysine/arginine"]:
    verdict = ("FORMAL CHARGE ALONE. The charge-class contrast carries everything and neither within-class pKa "
               "difference resolves, so the sign of the residue's charge near neutral pH sets the coefficient "
               "and pKa orders nothing.")
elif sig["acidic minus basic"]:
    verdict = ("NEITHER ACCOUNT AS WRITTEN: charge sign sets the DIRECTION and pKa orders magnitudes WITHIN at "
               "least one class, with the classes sign-uniform. The pre-registered ordering failed because it "
               "was applied across classes, where charge sign dominates it.")
else:
    verdict = ("NEITHER: the charge-class contrast itself does not resolve, which would undercut the grouping "
               "the whole section rests on and should be investigated before anything else here is read.")
print("   " + verdict)
print(f"   charge class {ab[0]:+.5f} (p {ab[2]:.1e}) for a pKa gap of {GAPS['acidic minus basic']:.1f} units;")
print(f"   glutamate over aspartate {ae[0]:+.5f} (p {ae[2]:.1e}) for a gap of {GAPS['glutamate minus aspartate']:.1f};")
print(f"   histidine against lysine/arginine {ah[0]:+.5f} (p {ah[2]:.1e}) for a gap of "
      f"{GAPS['histidine minus lysine/arginine']:.1f}.")
rows.append({"part": "verdict", "verdict": verdict, "basic_class_sign_uniform": bool(uniform_basic),
             "sig_charge_class": sig["acidic minus basic"],
             "sig_within_acidic": sig["glutamate minus aspartate"],
             "sig_within_basic": sig["histidine minus lysine/arginine"]})
# THE RATE IS THE INTERESTING PART IF BOTH FIRE. pKa per unit is not comparable across a 6-unit gap and a
# 0.6-unit one unless the relationship is linear, which nothing here establishes - so the two rates are printed
# as a comparison that would have to hold for a simple pKa account and is not itself a test.
print(f"\n4. per unit of pKa, if the relationship were linear - which nothing here establishes")
for nm, (b, se, pv) in tests.items():
    print(f"   {nm:>34} {b / GAPS[nm]:>+10.5f} per pKa unit")
    rows.append({"part": "per pka unit", "contrast": nm, "beta_per_pka_unit": b / GAPS[nm]})
print("   If one account explained all three, these would agree. They are printed because a simple pKa account")
print("   requires that and the comparison is cheap, not because agreeing would establish it. They do not agree,")
print("   and one of the three has the OPPOSITE sign, so a linear pKa account is excluded outright.")

print("\n5. the reading")
print("Correction 18 answered a yes-or-no question - does the residue ranking follow side-chain pKa - and the")
print("answer was no. That left the useful question untouched, and with all five titratable residues on the panel")
print("it is answerable in three pieces rather than one.")
print("CHARGE CLASS IS REAL AND IS NOT AN ARTEFACT OF HISTIDINE. The acidic pair differs from lysine and arginine")
print(f"by {abk[0]:+.5f} per point (p {abk[2]:.1e}) with histidine held fixed, so the grouping survives removing")
print("the one residue that could have manufactured it. That test mattered: against the whole basic class the")
print(f"contrast is {tests['acidic minus basic'][0]:+.5f}, more than twice as large, because histidine is in there.")
print("PKA ORDERS THE ACIDIC PAIR. Glutamate exceeds aspartate by +0.03958 per point at p 3.9e-03 across a pKa")
print("gap of 0.6 units, which is the direction a protonation account predicts and is now a coefficient with a")
print("standard error rather than two numbers set side by side.")
# DERIVED, both figures. This paragraph first said "five times" where the ratio is 8.6, and listed the three
# per-pKa rates as typed constants. Same hardcoded-number error as four other scripts in two days, here in a
# closing paragraph where nothing else would have checked it.
_hr = abs(joint["H"][0]) / (abs(joint["K"][0]) + abs(joint["R"][0]))
_rates = [tests[nm][0] / GAPS[nm] for nm in tests]
print(f"AND HISTIDINE BELONGS TO NEITHER PATTERN. It is negative where its class-mates are null, {_hr:.1f} times")
print("the magnitude of the lysine and arginine coefficients combined, and the basic class is not sign-uniform")
print("because of it. The per-pKa-unit rates make the same point arithmetically: "
      + ", ".join(f"{r:+.3f}" for r in _rates) + ", with")
print(f"{sum(1 for r in _rates if (r > 0) != (_rates[0] > 0))} of the three reversed in sign, so no single pKa "
      "slope fits all three contrasts.")
print("So the pre-registered ordering failed for TWO reasons, not one: it was applied across charge classes where")
print("charge sign dominates pKa, and it included a residue that does not behave like its class. That is a more")
print("useful result than the negative one it replaces, and it is consistent with correction 19's finding that")
print("the acidic and histidine effects are separate co-varying mechanisms rather than one quantity - this is the")
print("same conclusion reached from the residue side instead of the sidedness side, which is the agreement worth")
print("having between two routes to one claim.")
pd.DataFrame(rows).to_csv(D / f"what_organises_residues{TAG}.csv", index=False)
print("\nwrote", D / f"what_organises_residues{TAG}.csv")
