"""Can the recommended second forward arm actually be built? -> DATA_DIR/second_arm_editability.csv
docs/CHARGE_REVERSAL_DESIGN.md section 4b and docs/PREREGISTRATION.md section 6 both recommend adding AmtB (P69681, PF00909) as a second forward arm,
on the grounds that it is the most tractable of the three references that combine a high panel reliability with a positive coefficient - the cheapest
way to stop the whole experiment resting on one Pfam family's reference structures. That recommendation was made on PANEL-LEVEL merit and never checked
against the only thing that decides whether an arm can be built: whether the reference carries enough acidic residues at exposed periplasmic positions
to reverse.
It does not. Running the per-position scan for AmtB (scripts/37a, whose reference set is now configurable for exactly this purpose) gives 40 exposed
periplasmic positions of which only 5 carry D or E in the reference, and applying the design's own selection criterion from scripts/37b - a basic residue
tolerated at that position in at least 2% of bacteria - leaves 2. Both have a coefficient at or below zero, meaning the acidic state is associated with
a LOWER limit there, which is the opposite of the direction the arm would be built to test. Against the 14 positions the primary arm uses, an arm of 2
positions pointing the wrong way is not an arm.
The recommendation was a specific instance of the error scripts/54 was written to warn about: a reference's own coefficient correlates +0.817 with how
cleanly it can be scored, so panel-level merit largely measures annotation quality. But the tempting version of that criticism is too strong and is not
what the data say. Reliability and the count of reference D/E positions correlate +0.591 across the panel, so panel merit IS a moderate guide to
buildability and AmtB was a poor pick within a positive relationship rather than proof the two are unrelated. The error was narrower and more
embarrassing than 'used the wrong criterion': it was recommending an arm without running a filter that costs nothing, needs no compute, and would have
ranked AmtB seventh.
So this script does two things.
  1. THE VERDICT ON AMTB, from the full scan, against the design's own criteria rather than a looser reading of them.
  2. A CHEAP PRE-FILTER over all 32 membrane references that would have caught it before anything was recommended, and which now ranks the alternatives.
     Counting acidic residues at a reference's own exposed periplasmic positions needs nothing but the class table - no alignment, no scan, no compute -
     so there was never a reason not to run it first. Polymorphism still needs the scan, and that is stated: this filter establishes the ceiling on how
     many positions an arm could have, not how many would pass selection.
Usage: python3 scripts/71_second_arm_editability.py   (env DATA_DIR=results/broad2)
       the AmtB scan it reads is produced by:
       REFS='{"P69681": ["P69681", "PF00909"]}' CLASS_CSV=results/structure_panel_wide_classes.csv OUT_TAG=_amtb python3 scripts/37a_candidate_positions.py"""
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
BASIC_TOLERATED = 0.02      # the criterion scripts/37b applies: a basic residue present in >=2% of bacteria
PRIMARY_ARM_POSITIONS = 14  # what the NuoL + NuoN combined variant reverses, for scale
rows = []

print("1. the verdict on AmtB, from the full per-position scan")
scan = pd.read_csv(D / "candidate_positions_amtb.csv")
de = scan[scan.ref_aa.isin(["D", "E"])]
elig = de[de.frac_basic >= BASIC_TOLERATED]
print(f"   {len(scan)} exposed periplasmic positions; {len(de)} carry D or E in the reference; "
      f"{len(elig)} also tolerate a basic residue in >={BASIC_TOLERATED:.0%} of bacteria")
print(f"   {'pos':>5} {'aa':>3} {'n_aligned':>10} {'frac_basic':>11} {'beta':>8} {'p':>8}")
for r in de.sort_values("frac_basic", ascending=False).itertuples():
    mark = "  <- eligible" if r.frac_basic >= BASIC_TOLERATED else ""
    b = f"{r.beta:+.3f}" if np.isfinite(r.beta) else "n/a"
    p = f"{r.p:.3f}" if np.isfinite(r.p) else "n/a"
    print(f"   {int(r.pos):>5} {r.ref_aa:>3} {int(r.n_aligned):>10} {r.frac_basic:>11.3f} {b:>8} {p:>8}{mark}")
wrong_way = int((elig.beta <= 0).sum())
print(f"   of the {len(elig)} eligible, {wrong_way} have a coefficient at or below zero - the acidic state associated with a LOWER")
print(f"   limit, which is the opposite of the direction this arm would be built to test.")
print(f"   The primary arm reverses {PRIMARY_ARM_POSITIONS} positions. An arm of {len(elig)} is not an arm, and the recommendation is withdrawn.")
rows.append({"part": "amtb verdict", "ref": "P69681", "name": "AmtB", "peri_exp_positions": len(scan),
             "reference_DE": len(de), "eligible": len(elig), "eligible_wrong_sign": wrong_way,
             "primary_arm_positions": PRIMARY_ARM_POSITIONS})

print("\n2. the pre-filter that would have caught this, over all 32 membrane references")
print("   Counting acidic residues at a reference's own exposed periplasmic positions needs only the class table - no")
print("   alignment, no scan, no compute - so there was never a reason not to run it before recommending an arm.")
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
    out.append({"ref": str(r.label), "name": str(r.name)[:40], "peri_exp": len(aa),
                "reference_DE": int(np.isin(aa, ["D", "E"]).sum()),
                "reference_KR": int(np.isin(aa, ["K", "R"]).sum()),
                "reliability": float(rank.reliability.get(str(r.label), np.nan)),
                "beta_per_sd": float(rank.beta_per_sd.get(str(r.label), np.nan))})
F = pd.DataFrame(out).sort_values("reference_DE", ascending=False)
print(f"\n   {'ref':>8} {'name':40s} {'peri_exp':>9} {'D/E':>5} {'K/R':>5} {'reliability':>12} {'beta':>8}")
for x in F.head(12).itertuples():
    print(f"   {x.ref:>8} {x.name:40s} {x.peri_exp:>9} {x.reference_DE:>5} {x.reference_KR:>5} "
          f"{x.reliability:>12.3f} {x.beta_per_sd:>+8.3f}")
amtb_rank = int((F.reference_DE > F[F.ref == "P69681"].reference_DE.iloc[0]).sum()) + 1
print(f"\n   AmtB has {int(F[F.ref == 'P69681'].reference_DE.iloc[0])} reference D/E positions and ranks {amtb_rank} of {len(F)} on that count.")
_r = float(F[["reliability", "reference_DE"]].corr().iloc[0, 1])
print(f"   Reliability and D/E count correlate r = {_r:+.3f} across the panel. That is NOT the dismissal it would be")
print("   convenient to report: panel merit is a moderate guide to buildability, and AmtB was a poor pick within an overall")
print("   positive relationship rather than evidence the two properties are unrelated. The error was not using panel merit;")
print("   it was recommending an arm without running a filter that costs nothing and would have ranked AmtB 7th.")
rows += [{"part": "prefilter", **x} for x in F.to_dict("records")]
rows.append({"part": "prefilter summary", "amtb_DE_rank": amtb_rank, "n_refs": len(F),
             "r_reliability_vs_DE": float(F[["reliability", "reference_DE"]].corr().iloc[0, 1])})

print("\n3. what a replacement arm would need, and the candidates on this filter alone")
cand = F[(F.reference_DE >= 8) & (F.reliability >= 0.20)]
print(f"   Requiring at least 8 reference D/E positions (so that an arm is at least half the size of the primary one) and")
print(f"   reliability at or above the 0.20 floor scripts/54 established leaves {len(cand)} of {len(F)} references:")
if cand.empty:
    print("   NONE. On the panel as it stands there is no second forward arm that is both scorable and editable, which is a")
    print("   harder conclusion than withdrawing AmtB and should be recorded as such rather than softened.")
else:
    for x in cand.itertuples():
        print(f"     {x.ref} {x.name:40s} D/E {x.reference_DE:2d}  reliability {x.reliability:.3f}  beta {x.beta_per_sd:+.3f}")
    print("   Each still needs the per-position scan before any of them is recommended, because polymorphism is not in this")
    print("   filter. That scan is one command per reference and is the step that was skipped last time.")
rows.append({"part": "candidates", "n_passing": len(cand), "min_DE": 8, "min_reliability": 0.20,
             "refs": ";".join(cand.ref) if len(cand) else ""})

print("\n4. the replacement, scanned rather than recommended on the filter alone")
# the filter's top candidate, put through the same per-position scan AmtB was, which is the step skipped last time
cand_file = D / "candidate_positions_p24177.csv"
if cand_file.exists():
    p24 = pd.read_csv(cand_file)
    pde = p24[p24.ref_aa.isin(["D", "E"])]
    pel = pde[pde.frac_basic >= BASIC_TOLERATED]
    orig = pd.read_csv(D / "candidate_positions.csv")
    ode = orig[orig.ref_aa.isin(["D", "E"])]
    oel = ode[ode.frac_basic >= BASIC_TOLERATED]
    print(f"   P24177 (aminoglycoside efflux pump, PF00873): {len(p24)} exposed periplasmic positions, {len(pde)} carry D/E,")
    print(f"   {len(pel)} eligible - against {len(oel)} for YidC, NuoL and NuoN COMBINED ({len(ode)} D/E of {len(orig)} positions).")
    print(f"   Of the eligible, {int((pel.beta > 0).sum())} have a positive coefficient, against {int((oel.beta > 0).sum())} in the original set.")
    print(f"   {'pos':>5} {'aa':>3} {'frac_basic':>11} {'beta':>8} {'p':>7}")
    # stored as rows, not only printed: these per-position figures are quoted in the design document and the README, and
    # they live in the scan's output rather than this script's, so without this they are untraceable to the file that
    # cites them. scripts/69 caught exactly that.
    for r in pel.nlargest(5, "beta").itertuples():
        print(f"   {int(r.pos):>5} {r.ref_aa:>3} {r.frac_basic:>11.3f} {r.beta:>+8.3f} {r.p:>7.3f}")
        rows.append({"part": "replacement position", "ref": "P24177", "pos": int(r.pos), "ref_aa": r.ref_aa,
                     "frac_basic": float(r.frac_basic), "beta": float(r.beta), "p": float(r.p)})
    print(f"   So on buildability it is a better arm than the primary one, not merely an adequate second.")
    print(f"\n   One caution that belongs with it. The MEAN coefficient across its {len(pel)} eligible positions is "
          f"{pel.beta.mean():+.4f},")
    print("   essentially zero: the top positions are positive and the set as a whole is not, which is the same picture the")
    print("   per-position scan gave for the original panel (0 of 165 positions survive multiple-testing correction). Selecting")
    print("   positions by their own coefficient is therefore selecting largely on noise. That does not inflate the predicted")
    print("   effect, because the design sizes the arm from the COMPARTMENT-level composite rather than by summing per-position")
    print("   coefficients - but it does mean the choice of which positions to edit carries no evidential weight of its own, and")
    print("   an arm must be sized by the composite as it stands.")
    rows.append({"part": "replacement", "ref": "P24177", "name": "aminoglycoside efflux pump", "peri_exp_positions": len(p24),
                 "reference_DE": len(pde), "eligible": len(pel), "eligible_positive_beta": int((pel.beta > 0).sum()),
                 "mean_beta_eligible": float(pel.beta.mean()), "original_panel_eligible": len(oel),
                 "original_panel_positive": int((oel.beta > 0).sum())})
else:
    print(f"   {cand_file} not present; run the scan command in the docstring with OUT_TAG=_p24177 first")

pd.DataFrame(rows).to_csv(D / "second_arm_editability.csv", index=False)
print("\nwrote", D / "second_arm_editability.csv")
print("\nReading. The recommendation to add AmtB is withdrawn and P24177 replaces it. AmtB was chosen on panel-level reliability and")
print("coefficient, which is a moderate guide to buildability (r +0.591) rather than no guide at all - so the error was not the")
print("criterion but skipping a filter that costs nothing, needs no compute, and would have ranked AmtB seventh of 32. The filter")
print("belongs ahead of any future arm choice, and the per-position scan behind it belongs ahead of any recommendation.")
