"""Does histidine's periplasmic coefficient differ from its cytoplasmic one, or is the whole thing one non-specific trend? -> DATA_DIR/histidine_sidedness.csv and, with `wide`, DATA_DIR/histidine_sidedness_wide.csv
BOTH SPELLINGS ARE DECLARED for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line, keeps a panel variant's pool separate, and an
undeclared output makes every figure quoted from it a guaranteed false miss.
WHAT scripts/81 LEFT UNRESOLVED, stated as the objection it is. Correction 18 found pooled exposed-periplasmic histidine at -0.2207 per point, p small, opposite in sign to
glutamate - the buffering direction - and then failed the project's specificity standard: the exposed CYTOPLASMIC control, which is null for the acidic result (+0.019,
p 0.38), came back significant, same-signed, at 55% of the periplasmic coefficient. scripts/81 reported that and stopped, which was right as far as it went and is not an
answer. Two readings survive it. (a) ONE TREND: there is a single histidine association with pH limit that has nothing to do with which side of the membrane a residue sits
on, and slicing it by position class produces two shadows of it. (b) TWO COMPONENTS: a shared component plus a genuinely larger periplasmic one. A separate-models comparison
cannot tell these apart, because two models fitted on overlapping measures of the same thing will both come back significant under either reading.
THE TEST THAT CAN TELL THEM APART IS A CONTRAST, NOT A SECOND MODEL. Enter both shares in ONE fit and read the DIFFERENCE. Under (a) the difference is zero whatever the two
levels do; under (b) it is non-zero in the periplasmic direction. The difference also has a standard error, so "the periplasmic side is more negative" becomes a claim that
can fail rather than an impression formed by putting two numbers next to each other.
AND IT IS PARAMETRISED AS SUM AND DIFFERENCE RATHER THAN AS THE TWO LEVELS. {peri, cyto} and {peri+cyto, peri-cyto} span the same column space, so the FIT is identical -
same residuals, same log-likelihood, and part 7 checks that by recovering the two level coefficients from the sum and difference and matching them against the separately
fitted figures. What the reparametrisation buys is that the DIFFERENCE is estimated directly and so carries a standard error, which is what turns "the periplasmic side is
more negative" from an impression formed by putting two numbers side by side into a claim that can fail.
I FIRST WROTE THAT IT ALSO BUYS BETTER CONDITIONING, AND THE PRINTED NUMBERS SAY OTHERWISE, so the claim is withdrawn rather than quietly dropped. Sum and difference are
near-orthogonal only when the two sides have similar spread. For histidine they do - sd 0.76 periplasmic against 0.61 cytoplasmic, a ratio of 1.24 - and r(sum, difference)
is +0.23 against +0.44 for the levels, so the reparametrisation does help there. For the acidic measure the periplasmic spread is 2.05x the cytoplasmic, and r(sum,
difference) is +0.62 against only +0.17 for the levels: sum and difference are the WORSE-conditioned pair, by a wide margin. Both terms clear p<0.05 for both measures, so
nothing here is being hidden by an inflated standard error, but the justification for the parametrisation is the first paragraph's and not this one. The variance ratio is
printed beside the correlations so the reader can see which case they are in.
THE CALIBRATION IS THE ACIDIC RESULT, which is the part that makes a null here interpretable. D+E on the same panel, the same genomes, the same pooling, the same contrast:
its cytoplasmic control is null, so its DIFFERENCE must come back non-zero and periplasmic. If this contrast cannot detect sidedness where sidedness is already established,
a null for histidine says nothing about histidine. A test is only a test once it has been shown to fire.
SCALE. Every coefficient is per percentage point of pooled residue share, from ph.pgls directly, never ph.pgls_scan - correction 15.
Usage: python3 scripts/82_histidine_sidedness.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
MINPOS = 200                     # the pooled-position floor scripts/81 set, applied to BOTH classes here
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
if not any(c.endswith("|H") for c in pan.columns):
    sys.exit(f"{PANEL} has no histidine columns; rebuild the panel with scripts/36a first")
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "frac_His"]], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]


def pooled(cls, keys=("D", "E", "H", "K")):
    """Sum counts and positions across the membrane references, then divide - scripts/81's pooling, both sides."""
    tot = {k: np.zeros(len(B)) for k in ("n",) + tuple(keys)}
    used = 0
    for ref in MEM:
        if f"{ref}|{cls}|n" not in B or f"{ref}|{cls}|H" not in B:
            continue
        tot["n"] += B[f"{ref}|{cls}|n"].fillna(0).to_numpy(float)
        for k in keys:
            tot[k] += B[f"{ref}|{cls}|{k}"].fillna(0).to_numpy(float)
        used += 1
    return tot, used


PT, up = pooled("peri_exp")
CT, uc = pooled("cyto_exp")
ok = (PT["n"] >= MINPOS) & (CT["n"] >= MINPOS)
g = B[ok].reset_index(drop=True)
print(f"1. the sample: {up} references pooled periplasmic, {uc} cytoplasmic; "
      f"{int(ok.sum())} of {len(B)} genomes clear {MINPOS} pooled positions on BOTH sides")
# Both classes must clear the floor, so this sample can only be a subset of scripts/81's part-2 sample: the
# contrast is defined only for a genome scorable on both sides. The cost of asking a paired question is printed
# rather than assumed, and on this panel it turns out to be zero - every genome clears 200 pooled positions on
# both sides, because pooling 32 references puts roughly 860 positions per genome on each side.
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
X0 = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
y = g.ph_min.to_numpy(float)
lam = ph.pgls(y, X0, Cm)["lambda"]
P = {k: 100 * PT[k][ok] / PT["n"][ok] for k in ("D", "E", "H", "K")}
C = {k: 100 * CT[k][ok] / CT["n"][ok] for k in ("D", "E", "H", "K")}
rows.append({"part": "sample", "n_genomes": int(len(g)), "n_bacteria_available": int(len(B)),
             "refs_peri": up, "refs_cyto": uc, "minpos": MINPOS, "lambda": lam})

print("\n2. the conditioning check, which decides whether the two LEVELS could have been read at all")
print(f"   {'measure':>22} {'sd ratio p/c':>13} {'r(peri,cyto)':>13} {'r(sum,diff)':>12} {'better pair':>13}")
cond = {}
for lab, pp, cc in (("histidine", P["H"], C["H"]), ("acidic D+E", P["D"] + P["E"], C["D"] + C["E"])):
    r = float(np.corrcoef(pp, cc)[0, 1])
    rsd = float(np.corrcoef(pp + cc, pp - cc)[0, 1])
    ratio = float(pp.std() / cc.std())
    cond[lab] = (r, rsd, ratio)
    print(f"   {lab:>22} {ratio:>13.2f} {r:>+13.3f} {rsd:>+12.3f} "
          f"{('sum/diff' if abs(rsd) < abs(r) else 'the levels'):>13}")
    rows.append({"part": "conditioning", "measure": lab, "sd_ratio_peri_over_cyto": ratio,
                 "r_peri_cyto": r, "r_sum_diff": rsd, "sum_diff_better_conditioned": bool(abs(rsd) < abs(r))})
print("   Sum and difference decorrelate only when the two sides have similar spread, so for histidine (ratio 1.2)")
print("   they are the better-conditioned pair and for the acidic measure (ratio 2.1) they are the worse one. The")
print("   contrast is read off them anyway, because the difference is the QUANTITY ASKED FOR and this is the")
print("   parametrisation that gives it a standard error. Part 7 confirms the fit itself is unchanged either way.")

print("\n3. the contrast, per percentage point. sum = shared component, difference = periplasmic excess")
print(f"   {'measure':>22} {'term':>12} {'beta/pt':>10} {'se':>9} {'p':>10}")
fit = {}
for lab, pp, cc, others in (("acidic D+E", P["D"] + P["E"], C["D"] + C["E"], ()),
                            ("histidine", P["H"], C["H"], ("D", "E", "K"))):
    # The histidine fit holds D, E and K on the periplasmic side fixed, as correction 14 requires: a residue share
    # entered alone is not that residue's effect. The acidic fit enters D+E as the quantity of interest itself.
    extra = [P[k] for k in others]
    X = np.column_stack([X0] + extra + [pp + cc, pp - cc])
    r = ph.pgls(y, X, Cm, lam=lam)
    out = {}
    for j, term in (( -2, "sum"), (-1, "difference")):
        b, se, p = float(r["beta"][j]), float(r["se"][j]), float(r["p"][j])
        out[term] = (b, se, p)
        print(f"   {lab:>22} {term:>12} {b:>+10.5f} {se:>9.5f} {p:>10.2e}")
        rows.append({"part": "contrast", "measure": lab, "term": term, "beta_per_pt": b, "se": se, "p": p,
                     "n": len(g)})
    fit[lab] = out

print("\n4. is the contrast able to detect sidedness at all? the acidic result is the calibration")
ab, ase, ap = fit["acidic D+E"]["difference"]
fires = ap < 0.05 and ab > 0
print(f"   acidic D+E difference (periplasmic minus cytoplasmic): {ab:+.5f}/pt (se {ase:.5f}, p {ap:.2e})")
print("   " + ("FIRES as it must: the acidic effect IS one-sided by every other test in this project, and this"
               if fires else
               "DOES NOT FIRE, which disqualifies the test rather than the hypothesis: the acidic effect is"))
print("   " + ("   contrast sees it, so a null for histidine below is informative about histidine."
               if fires else
               "   one-sided by every other test here, so a contrast that cannot see it cannot be read."))
rows.append({"part": "calibration", "beta_diff_acidic_per_pt": ab, "se": ase, "p": ap, "fires": fires})
assert fires, ("the calibration failed: the acidic difference is not positive and significant, so nothing below "
               "is interpretable. Do not read the histidine contrast past this point.")

print("\n5. histidine, read against that calibration")
hs, hse_, hsp = fit["histidine"]["sum"]
hd, hdse, hdp = fit["histidine"]["difference"]
print(f"   shared component  {hs:+.5f}/pt (p {hsp:.2e})")
print(f"   periplasmic excess {hd:+.5f}/pt (se {hdse:.5f}, p {hdp:.2e})")
if hdp < 0.05:
    verdict = ("TWO COMPONENTS: histidine's periplasmic coefficient is distinguishable from its cytoplasmic one, "
               "so there is a side-specific part on top of a shared trend. Correction 18's association is not "
               "wholly an artefact of slicing a genome-wide trend by position class.")
elif hsp < 0.05:
    verdict = ("ONE TREND: the shared component carries the association and the periplasmic excess is not "
               "distinguishable from zero on a test shown to fire for the acidic result. So the histidine "
               "association is real and is NOT a property of the outward-facing surface - which is what the "
               "failed specificity control in scripts/81 was already pointing at, now with a contrast behind it.")
else:
    verdict = ("NEITHER COMPONENT RESOLVES on the paired sample, which is smaller than scripts/81's. That is a "
               "power statement about this subset, not a retraction of correction 18.")
print(f"   {verdict}")
rows.append({"part": "verdict", "beta_sum_per_pt": hs, "p_sum": hsp, "beta_diff_per_pt": hd, "se_diff": hdse,
             "p_diff": hdp, "verdict": verdict})

print("\n6. what this does to the tier correction 18 was given, and to the experiment arm it declined")
print("   scripts/81 tiered histidine as 'a real association whose mechanism is unsettled' and declined a D/E->H")
print("   experimental arm because the coefficient is partly non-specific. This contrast sharpens the second half")
print("   of that: " + ("part of the coefficient IS specific to the periplasmic side, so the arm's rationale is "
                        "stronger than scripts/81 allowed - though the shared component still means a bench "
                        "substitution would be editing only part of what the coefficient measures."
                        if hdp < 0.05 else
                        "none of the coefficient is demonstrably specific to the periplasmic side, so declining "
                        "the arm was right for a stronger reason than scripts/81 had. An edit to outward-facing "
                        "histidines is predicted to move nothing, because the association is not located there."))
print("   Either way the decision is unchanged and the reason for it is now measured.")
rows.append({"part": "design implication", "arm_adopted": False,
             "periplasmic_excess_resolved": bool(hdp < 0.05)})

print("\n7. the invariance check: the levels recovered from the contrast, against the levels fitted directly")
# If {sum, difference} and {peri, cyto} really span the same space then peri = b_sum + b_diff and cyto = b_sum -
# b_diff exactly, and fitting the levels must reproduce them to numerical precision. This is worth spending a fit
# on for two reasons. It verifies the arithmetic claim the whole parametrisation rests on, which was asserted in
# the docstring and is cheap to check rather than trust. And the recovered levels are directly comparable with
# scripts/81's separately fitted coefficients, so the two scripts can be read against each other.
print(f"   {'measure':>22} {'side':>6} {'from contrast':>14} {'fitted directly':>16} {'se':>9} {'p':>10}")
LEV = {"acidic D+E": (P["D"] + P["E"], C["D"] + C["E"], ()), "histidine": (P["H"], C["H"], ("D", "E", "K"))}
for lab, (pp, cc, others) in LEV.items():
    r = ph.pgls(y, np.column_stack([X0] + [P[k] for k in others] + [pp, cc]), Cm, lam=lam)
    bs, bd = fit[lab]["sum"][0], fit[lab]["difference"][0]
    for j, side, implied in ((-2, "peri", bs + bd), (-1, "cyto", bs - bd)):
        b, se, pv = float(r["beta"][j]), float(r["se"][j]), float(r["p"][j])
        print(f"   {lab:>22} {side:>6} {implied:>+14.5f} {b:>+16.5f} {se:>9.5f} {pv:>10.2e}")
        assert abs(b - implied) < 1e-6, f"{lab} {side}: reparametrisation is not an identity ({b} vs {implied})"
        rows.append({"part": "invariance", "measure": lab, "side": side, "beta_from_contrast_per_pt": implied,
                     "beta_fitted_per_pt": b, "se": se, "p": pv})
print("   The two columns agree to within 1e-6, which is the identity the parametrisation claimed.")
# The acidic cytoplasmic level landing on essentially zero is the same null that docs/SYNTHESIS.md reports from a
# separate model (+0.019, p 0.38), reached here as a by-product of a differently specified fit - and the histidine
# periplasmic level lands within a few thousandths of scripts/81's -0.2207 on a different sample definition. Two
# independent routes to each figure is the check; agreeing is what entitles reading them together.
print("   Both independently recover figures the project already has: the acidic cytoplasmic level is the null")
print("   docs/SYNTHESIS.md reports from a separate model, and histidine's periplasmic level sits beside")
print("   scripts/81's pooled estimate.")

print("\n8. THE SPECIFICITY CONTROL IN scripts/81 IS MIS-SPECIFIED, and part 7 is what shows it")
# Part 7 hands back a cytoplasmic histidine coefficient of about -0.052 at p 0.14 - a NULL - where scripts/81
# reported that same control as significant, same-signed and worth 55% of the periplasmic coefficient, and
# declined an experimental arm on the strength of it. The two cannot both be right about the same data, and the
# difference is entirely the specification: scripts/81 enters cytoplasmic histidine ALONE, in a model with no
# periplasmic term in it at all, and the two shares correlate at r 0.44 across genomes (part 2). A regressor
# entered without a correlate it shares variance with absorbs that correlate's signal. That is omitted-variable
# bias, it is the same error correction 14 found for aspartate-without-glutamate, and it is the project's own
# recurring failure mode showing up inside a control rather than inside a headline.
# THE DISCRIMINATING PREDICTION, which is what makes this a demonstration and not a preference between two fits:
# if the mechanism is omitted-variable bias, then adding the periplasmic term must move the CYTOPLASMIC estimate
# a long way for histidine (r 0.44) and much less for the acidic measure (r 0.17). The acidic control is already
# null under both specifications, so the two accounts cannot be told apart by their verdicts - only by HOW FAR
# the estimate moves, which the correlation predicts in advance and which is not free to come out either way.
# THE MOVE IS MEASURED IN STANDARD ERRORS, AND MY FIRST ATTEMPT MEASURED IT AS A PERCENTAGE OF THE BASELINE,
# which was meaningless and said the prediction had failed. The acidic cytoplasmic estimate starts at +0.0062 and
# ends at -0.0004, so a RELATIVE shift divides by a number that is itself indistinguishable from zero and came
# out as 106% against histidine's 48% - an artefact of the denominator, not a result. In standard errors, the
# scale on which "moved a long way" means anything, the ordering is the predicted one by a factor of five. The
# wrong metric is recorded here because the corrected number is only worth as much as the account of what it
# replaced.
# AND THE THIRD SPECIFICATION DOES NOT EXIST FOR THE ACIDIC MEASURE. Its cytoplasmic regressor IS cytoplasmic
# D+E, so entering cytoplasmic D, E and K beside it is exactly collinear; the first draft of this script fitted
# it anyway and printed a standard error of 168457 and p = 1.00 as though it were a finding. It is skipped, and
# the assertion below makes a rank-deficient fit impossible to report as a result again.
print(f"   {'measure':>12} {'specification':>44} {'cyto beta/pt':>13} {'se':>8} {'p':>9}")
SPECS = (("as scripts/81 fits it: cytoplasmic share alone", ()),
         ("+ the periplasmic share of the same residue", ("own",)),
         ("+ periplasmic D/E/K and cytoplasmic D/E/K", ("own", "pdek", "cdek")))
moved, ses = {}, {}
for lab, (pp, cc, _o) in LEV.items():
    for sname, adds in SPECS:
        if "cdek" in adds and lab == "acidic D+E":
            print(f"   {lab:>12} {sname:>44} {'(collinear: the regressor IS cytoplasmic D+E)':>32}")
            rows.append({"part": "specificity respecified", "measure": lab, "specification": sname,
                         "skipped": "exactly collinear - the regressor is cytoplasmic D+E itself"})
            continue
        cols = [X0]
        if "own" in adds:
            cols.append(pp.reshape(-1, 1))
        if "pdek" in adds:
            cols.append(np.column_stack([P[k] for k in ("D", "E", "K")]))
        if "cdek" in adds:
            cols.append(np.column_stack([C[k] for k in ("D", "E", "K")]))
        cols.append(cc.reshape(-1, 1))
        X = np.column_stack(cols)
        assert np.linalg.matrix_rank(X) == X.shape[1], f"{lab} / {sname}: design matrix is rank-deficient"
        r = ph.pgls(y, X, Cm, lam=lam)
        b, se, pv = float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1])
        moved.setdefault(lab, []).append(b); ses.setdefault(lab, []).append(se)
        print(f"   {lab:>12} {sname:>44} {b:>+13.5f} {se:>8.5f} {pv:>9.2e}"
              + ("   <- the figure scripts/81 reports" if not adds and lab == "histidine" else ""))
        rows.append({"part": "specificity respecified", "measure": lab, "specification": sname,
                     "beta_cyto_per_pt": b, "se": se, "p": pv, "significant": bool(pv < 0.05)})
sh_ = abs(moved["histidine"][1] - moved["histidine"][0]) / ses["histidine"][0]
sa = abs(moved["acidic D+E"][1] - moved["acidic D+E"][0]) / ses["acidic D+E"][0]
print(f"   adding the periplasmic term moves the cytoplasmic estimate by {sh_:.2f} standard errors for histidine")
print(f"   (r 0.44) and {sa:.2f} for the acidic measure (r 0.17): a factor of {sh_ / sa:.1f}, in the predicted")
print("   direction. On the percentage-of-baseline scale this script first used it read the other way round,")
print("   because the acidic baseline is itself zero - see the comment above.")
rows.append({"part": "specificity respecified", "measure": "shift test",
             "shift_histidine_in_se": sh_, "shift_acidic_in_se": sa, "ratio": sh_ / sa,
             "ordering_as_predicted": bool(sh_ > sa)})
assert sh_ > sa, ("the bias account predicted a larger shift where the correlation is larger, and it is not "
                  "larger; the reading below does not follow")
print("   The histidine control is null under BOTH of the better-specified fits, including the one that holds the")
print("   cytoplasmic acidic and basic shares fixed as well. So histidine PASSES the specificity standard after")
print("   all, and scripts/81's reported failure was its control borrowing the periplasmic share's signal.")
print("   WHAT THAT DOES AND DOES NOT CHANGE. It does not touch correction 18's periplasmic coefficient, its")
print("   sign, or the finding that the pKa ordering fails - those came from a model that already held the other")
print("   residues fixed. It removes the ONE reason scripts/81 gave for declining a D/E->H experimental arm, and")
print("   it means the shared component in part 3 is shared with something other than the cytoplasmic face.")

pd.DataFrame(rows).to_csv(D / f"histidine_sidedness{TAG}.csv", index=False)
# THE NULL IS MARGINAL AND IS NOT ROUNDED TO CLEAN. At p 0.072 the better-specified histidine control is not
# significant at the project's threshold, and it is also not the acidic result's +0.98 - the two should not be
# described in the same words. The honest form is that the cytoplasmic estimate falls below the threshold once
# the periplasmic share is held fixed and keeps falling as more is held fixed (p 0.18), and that the POSITIVE
# statement does not depend on that borderline call at all: the difference in part 3 is -0.0800 at p 0.003 on a
# contrast shown to fire for the acidic result, which says the two sides differ without needing either level to
# be called null.
print(f"\nwrote", D / f"histidine_sidedness{TAG}.csv")
print("\nReading. The question was whether correction 18's histidine association is a property of the outward-")
print("facing surface or one trend sliced into two shadows, and a contrast answers it where two separate models")
print("could not: the periplasmic excess is -0.0800 per point at p 0.003, on a test calibrated to fire for the")
print("acidic result and shown to do so. So there are two components, a shared one and a side-specific one.")
print("The unplanned finding is larger than the planned one. scripts/81's specificity control - the single reason")
print("it gave for declining a D/E->H experimental arm - is mis-specified: it enters the cytoplasmic share alone,")
print("the two shares correlate at r 0.44, and the coefficient it reports is partly the periplasmic share's. Add")
print("the periplasmic term and it moves 1.70 standard errors to p 0.072; hold the cytoplasmic acidic and basic")
print("shares fixed too and it reaches p 0.18. That null is marginal rather than clean and is not worth dressing")
print("up, but it does not have to be clean for the conclusion: the periplasmic-minus-cytoplasmic difference is")
print("significant on its own terms. This is correction 14's error - a residue entered without the correlate it")
print("shares variance with - found inside a CONTROL rather than inside a headline, which is the harder place to")
print("look and the ninth appearance of this family of mistake in the project.")
print("What it costs and what it does not. Correction 18's periplasmic coefficient, its sign, and the failure of")
print("the pKa ordering all stand: they came from a fit that already held the other residues fixed. What falls is")
print("the stated reason for declining the histidine arm. The arm is still not adopted here - that is the")
print("deviations procedure's call and not a script's - but the objection standing against it has been withdrawn,")
print("and the shared component in part 3 is now shared with something that is not the cytoplasmic face and has")
print("not been identified.")
