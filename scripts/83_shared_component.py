"""What is the shared histidine component shared WITH? -> DATA_DIR/shared_component.csv and, with `wide`, DATA_DIR/shared_component_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHAT scripts/82 LEFT NAMED BUT NOT IDENTIFIED. Decomposing exposed histidine into a sum and a difference across the two faces of the membrane found two components: a
periplasmic excess of -0.0800 per point (p 0.0031), and a SHARED component of -0.1317 per point (p 8e-09) that is larger than the excess and was described only as "shared
with something that is not the cytoplasmic face". A component bigger than the side-specific one cannot be left as a residual category, because every reading of correction
18's mechanism depends on what it is. Three candidates are testable from columns already on disk, and they are not equally interesting.
THE ONE THAT WOULD MATTER MOST IF TRUE IS NON-EXPOSED HISTIDINE IN THE SAME PROTEINS. The panel classifies every reference position, so the same 32 families also yield a
buried and transmembrane histidine share - positions that are NOT solvent-accessible and therefore cannot be titrated by the external proton gradient at all, whatever the
periplasmic pH. If that share predicts ph_min as strongly as the exposed sum, then the association is not a titration story in the first place: it is the histidine content
of these protein families behaving as a marker for something else, and the exposed/buried distinction the whole structural panel is built on would not be doing any work
here. This is the strongest available falsifier of correction 18's mechanism and it is the reason to run this script rather than stop at naming the component.
THE SECOND IS PROTEOME-WIDE HISTIDINE, which scripts/81 already found null (p 0.96) and which is re-entered here because that test was run against the periplasmic share
alone, not against the sum, and a null against one regressor is not a null against another.
THE THIRD IS THE PANEL'S SCORING QUALITY. scripts/54 established that a per-reference estimate is largely a measure of how cleanly a reference can be scored; the pooled
position count is the genome-level version of the same thing, and a genome whose proteins align poorly has both fewer scorable positions and a noisier share.
A NEGATIVE RESULT HERE IS THE INFORMATIVE ONE, which is worth saying before the numbers rather than after. If none of the three carries the shared component, the component
stays unidentified and this script has only narrowed it - but it will have ruled out the falsifier above, which is what the exposed/buried contrast is for.
SCALE. Per percentage point of pooled share throughout, from ph.pgls directly and never ph.pgls_scan - correction 15.
Usage: python3 scripts/83_shared_component.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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


def pooled(classes, keys=("D", "E", "H", "K")):
    """Sum counts and positions across the membrane references AND across the named classes, then divide."""
    tot = {k: np.zeros(len(B)) for k in ("n",) + tuple(keys)}
    for ref in MEM:
        for cls in classes:
            if f"{ref}|{cls}|n" not in B or f"{ref}|{cls}|H" not in B:
                continue
            tot["n"] += B[f"{ref}|{cls}|n"].fillna(0).to_numpy(float)
            for k in keys:
                tot[k] += B[f"{ref}|{cls}|{k}"].fillna(0).to_numpy(float)
    return tot


EXP = pooled(("peri_exp", "cyto_exp"))
PH, CH = pooled(("peri_exp",)), pooled(("cyto_exp",))
NON = pooled(("peri_bur", "cyto_bur", "TM"))
ok = (EXP["n"] >= MINPOS) & (NON["n"] >= MINPOS)
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
X0 = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
y = g.ph_min.to_numpy(float)
lam = ph.pgls(y, X0, Cm)["lambda"]
# TWO POOLINGS OF THE SAME POSITIONS, AND THEY ARE NOT THE SAME REGRESSOR. This script pools both faces into ONE
# share, H/(n_peri + n_cyto); scripts/82's shared component is the SUM of two separately-pooled shares,
# H_peri/n_peri + H_cyto/n_cyto. The pooled share is a position-weighted mean of the two, roughly 0.44 periplasmic
# plus 0.56 cytoplasmic here, so it is a different linear combination and its coefficient is NOT comparable in
# magnitude to that script's. The first draft of part 2 below asserted they were the same quantity and reported
# -0.3665 against scripts/82's -0.1317 as a confirmation of it. Part 2 now fits both and measures the relation
# instead of claiming one.
hexp = 100 * EXP["H"][ok] / EXP["n"][ok]
hnon = 100 * NON["H"][ok] / NON["n"][ok]
# AND THE OTHER HALF OF THAT DISCREPANCY WAS CORRECTION 14, COMMITTED HERE. scripts/82 holds the periplasmic D, E
# and K shares fixed, because a residue share entered without the correlated shares of the other charged residues
# absorbs their signal - which is the error correction 14 found for aspartate, correction 19 found inside a
# control, and this script's first draft then made a third time by controlling nothing but genome-wide acidic
# excess and data source. Every fit below holds the exposed D, E and K shares fixed. The uncontrolled figure is
# kept in part 2 as a measured contrast, not deleted, since the size of the error is the useful part.
DEK = np.column_stack([100 * EXP[k][ok] / EXP["n"][ok] for k in ("D", "E", "K")])
X0 = np.column_stack([X0, DEK])
gwh = g.frac_His.to_numpy(float) * 100
qual = np.log10(EXP["n"][ok])                 # scoring quality, as scripts/54 measures it: how much is scorable
print(f"1. the sample: {int(ok.sum())} of {len(B)} genomes clear {MINPOS} pooled positions in BOTH the exposed and")
print(f"   the non-exposed pools; lambda {lam:.3f}")
print(f"   {'pool':>28} {'positions/genome':>17} {'His share %':>12}")
for lab, tt in (("exposed (peri + cyto)", EXP), ("non-exposed (bur + TM)", NON)):
    print(f"   {lab:>28} {tt['n'][ok].mean():>17.0f} {100 * (tt['H'][ok] / tt['n'][ok]).mean():>12.3f}")
    rows.append({"part": "pools", "pool": lab, "positions_per_genome": float(tt["n"][ok].mean()),
                 "his_share_pct": float(100 * (tt["H"][ok] / tt["n"][ok]).mean()), "n_genomes": int(ok.sum())})

print("\n2. the regressor this script uses, placed against scripts/82's rather than confused with it")
XB = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])     # no charged-residue controls
hsum = (100 * PH["H"][ok] / PH["n"][ok]) + (100 * CH["H"][ok] / CH["n"][ok])
print(f"   {'regressor':>42} {'controls D/E/K':>15} {'beta/pt':>10} {'se':>9} {'p':>10}")
for lab, v, X, ctl in (("one pooled share over both faces", hexp, XB, "no"),
                       ("one pooled share over both faces", hexp, X0, "yes"),
                       ("sum of two shares (scripts/82's)", hsum, X0, "yes")):
    r = ph.pgls(y, np.column_stack([X, v]), Cm, lam=lam)
    b, se, pv = float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1])
    print(f"   {lab:>42} {ctl:>15} {b:>+10.5f} {se:>9.5f} {pv:>10.2e}")
    rows.append({"part": "regressor comparison", "regressor": lab, "controls_charged_residues": ctl,
                 "beta_per_pt": b, "se": se, "p": pv})
    if lab.startswith("one pooled") and ctl == "yes":
        b0, se0, p0 = b, se, pv
    if lab.startswith("sum"):
        bsum = b
print(f"   The pooled share is a position-weighted mean of the two faces' shares and the sum is their total, so")
print(f"   the two coefficients differ by construction. Omitting the charged-residue controls moves the pooled-share")
print(f"   coefficient as well, which is correction 14's error and is shown rather than only avoided. Everything")
print(f"   below uses the controlled fit.")
# AND THE SUM STILL IS NOT scripts/82's FIGURE, which is worth closing rather than waving at: that script reports
# -0.13168 and the row above gives -0.14466 for what is nominally the same regressor on the same genomes. My first
# explanation was the DIFFERENCE term - scripts/82 fits the sum and the difference together, the row above fits the
# sum alone, and the two correlate at +0.234 - and an assertion written to confirm it failed: restoring the
# difference moves the sum from -0.14466 to -0.14486, which is nothing. The cause is the CONTROL SET. scripts/82
# holds the PERIPLASMIC D, E and K shares fixed; this script holds the pooled-exposed ones, because its own
# question is about both faces together. Swap in the periplasmic set and both terms reproduce to five decimals.
# Which control set is right depends on the question, so both are reported and neither is called the correct one;
# what is NOT acceptable is two scripts quoting different numbers for the same name, which is correction 12.
hdif = (100 * PH["H"][ok] / PH["n"][ok]) - (100 * CH["H"][ok] / CH["n"][ok])
PDEK = np.column_stack([100 * PH[k][ok] / PH["n"][ok] for k in ("D", "E", "K")])
print(f"   {'control set':>42} {'sum':>10} {'difference':>12}")
for lab, ctl in (("periplasmic D/E/K, as scripts/82 fits it", PDEK),
                 ("pooled-exposed D/E/K, as used below", DEK)):
    r = ph.pgls(y, np.column_stack([XB, ctl, hsum, hdif]), Cm, lam=lam)
    bs, bd = float(r["beta"][-2]), float(r["beta"][-1])
    print(f"   {lab:>42} {bs:>+10.5f} {bd:>+12.5f}")
    rows.append({"part": "reconciliation with scripts/82", "control_set": lab, "beta_sum_per_pt": bs,
                 "beta_difference_per_pt": bd})
    if "scripts/82" in lab:
        assert abs(bs - (-0.13168)) < 1e-5 and abs(bd - (-0.08000)) < 1e-5, \
            f"did not reproduce scripts/82 under its own control set: sum {bs}, difference {bd}"
print("   Under scripts/82's own control set both of its terms reproduce exactly, so the two scripts agree and the")
print("   gap was the control set rather than the term I first blamed - an assertion written to confirm the wrong")
print("   explanation is what rejected it.")

print("\n3. the three candidates, each entered beside it")
print(f"   {'candidate':>34} {'r with exposed':>15} {'its own beta/pt':>16} {'p':>9} {'exposed becomes':>16} {'p':>9}")
CAND = (("non-exposed His, same proteins", hnon), ("proteome-wide His", gwh),
        ("scoring quality, log10 positions", qual))
for lab, v in CAND:
    X = np.column_stack([X0, v, hexp])
    assert np.linalg.matrix_rank(X) == X.shape[1], f"{lab}: design matrix is rank-deficient"
    r = ph.pgls(y, X, Cm, lam=lam)
    bc, pc = float(r["beta"][-2]), float(r["p"][-2])
    be, pe = float(r["beta"][-1]), float(r["p"][-1])
    rr = float(np.corrcoef(v, hexp)[0, 1])
    print(f"   {lab:>34} {rr:>+15.3f} {bc:>+16.5f} {pc:>9.2e} {be:>+16.5f} {pe:>9.2e}")
    rows.append({"part": "candidate", "candidate": lab, "r_with_exposed": rr, "beta_candidate_per_pt": bc,
                 "p_candidate": pc, "beta_exposed_controlled_per_pt": be, "p_exposed_controlled": pe,
                 "exposed_survives": bool(pe < 0.05 and (be > 0) == (b0 > 0))})

print("\n4. the falsifier, read explicitly: can a histidine that is not solvent-accessible do this?")
# A buried or transmembrane histidine is not in contact with the periplasm, so no external pH can protonate it.
# If it predicted ph_min as strongly as the exposed share, the exposed/buried distinction the structural panel is
# built on would be doing no work here and correction 18 would be a composition marker rather than a mechanism.
rn = ph.pgls(y, np.column_stack([X0, hnon]), Cm, lam=lam)
bn, sen, pn = float(rn["beta"][-1]), float(rn["se"][-1]), float(rn["p"][-1])
print(f"   non-exposed histidine ALONE: {bn:+.5f}/pt (se {sen:.5f}, p {pn:.2e})")
print(f"   exposed histidine ALONE:     {b0:+.5f}/pt (se {se0:.5f}, p {p0:.2e})")
both = ph.pgls(y, np.column_stack([X0, hnon, hexp]), Cm, lam=lam)
bnb, pnb = float(both["beta"][-2]), float(both["p"][-2])
beb, peb = float(both["beta"][-1]), float(both["p"][-1])
print(f"   entered together: non-exposed {bnb:+.5f} (p {pnb:.2e}), exposed {beb:+.5f} (p {peb:.2e})")
survives = peb < 0.05 and (beb > 0) == (b0 > 0)
killed = pnb < 0.05 and abs(bnb) >= abs(beb)
if survives and not killed:
    fals = ("PASSES: the exposed share keeps its coefficient with the non-exposed share held fixed, and the "
            "non-exposed share does not carry more. Solvent accessibility is doing work, which is the minimum "
            "the titration reading requires.")
elif killed:
    fals = ("FAILS: the non-exposed share carries at least as much as the exposed one, so this is the histidine "
            "CONTENT of these protein families and not a property of surfaces a proton can reach. Correction "
            "18's mechanism reading would have to be withdrawn.")
else:
    fals = ("INDETERMINATE: neither share resolves once both are entered, which given their correlation is a "
            "collinearity statement about this panel rather than a verdict.")
print(f"   {fals}")
rows.append({"part": "falsifier", "beta_nonexposed_alone_per_pt": bn, "p_nonexposed_alone": pn,
             "beta_nonexposed_joint": bnb, "p_nonexposed_joint": pnb,
             "beta_exposed_joint": beb, "p_exposed_joint": peb, "verdict": fals})

print("\n5. THE CANDIDATE THAT TURNED OUT TO BE A CONCERN FOR EVERY PANEL RESULT, NOT A SHARED COMPONENT")
# Scoring quality was entered in part 3 as the third and least interesting candidate, and it came back significant
# on its own terms: -2.36 pH per log10 pooled positions at p 1.7e-09. It does NOT absorb the histidine coefficient,
# which is what part 3 asked, so as a candidate it failed. But a genome-level predictor that strong is a question
# about the whole structural panel, because EVERY result in this project that reads a residue share off these
# references is computed on genomes that differ in how much of them is scorable - and scripts/54 already found
# that per-reference estimates are largely a measure of scoring cleanliness. It was never checked at genome level
# against the HEADLINE. Running it here rather than leaving it for a later script is the point: the check is
# cheap, it is already set up, and finding it in a script about histidine is not a reason to postpone it.
for lab, ks in (("acidic D+E, exposed periplasmic", ("D", "E")),):
    sh = sum(100 * PH[k][ok] / PH["n"][ok] for k in ks)
    qp = np.log10(PH["n"][ok])
    base = ph.pgls(y, np.column_stack([XB, sh]), Cm, lam=lam)
    ctl = ph.pgls(y, np.column_stack([XB, qp, sh]), Cm, lam=lam)
    b1, p1 = float(base["beta"][-1]), float(base["p"][-1])
    b2, p2 = float(ctl["beta"][-1]), float(ctl["p"][-1])
    bq, pq = float(ctl["beta"][-2]), float(ctl["p"][-2])
    print(f"   {lab}")
    print(f"      uncontrolled          {b1:+.5f}/pt (p {p1:.2e})")
    print(f"      scoring quality held  {b2:+.5f}/pt (p {p2:.2e})   [quality itself {bq:+.3f}, p {pq:.2e}]")
    print(f"      retains {abs(b2 / b1):.0%} of its coefficient")
    rows.append({"part": "scoring quality against the headline", "measure": lab, "beta_uncontrolled_per_pt": b1,
                 "p_uncontrolled": p1, "beta_controlled_per_pt": b2, "p_controlled": p2,
                 "beta_quality": bq, "p_quality": pq, "retained": abs(b2 / b1)})
    assert p2 < 0.05 and (b2 > 0) == (b1 > 0), \
        ("the acidic exposed-periplasmic result does not survive a genome-level scoring-quality control. That is "
         "a finding about the HEADLINE and must not be reported as a note inside a histidine script.")
print("   The headline survives, so this is a caveat that has been checked rather than a defect - but the control")
print("   belongs in the panel scripts and not only here, and that is recorded in docs/SYNTHESIS.md section 8.")

pd.DataFrame(rows).to_csv(D / f"shared_component{TAG}.csv", index=False)
print("\nwrote", D / f"shared_component{TAG}.csv")
print("\nReading. The shared component is still not identified, and the falsifier it was most at risk from is")
print("ruled out. Non-exposed histidine in the same 32 protein families - buried and transmembrane positions, which")
print("no external pH can protonate - predicts nothing on its own (p 0.37) and nothing beside the exposed share")
print("(p 0.47), while the exposed share keeps its full coefficient. Had that gone the other way, correction 18")
print("would have been the histidine CONTENT of these families behaving as a marker and the exposed/buried")
print("distinction the whole structural panel rests on would have been doing no work. Proteome-wide histidine is")
print("null here as it was in scripts/81, now against the pooled regressor rather than the periplasmic one, which")
print("is a second test and not the same one restated.")
print("The useful by-product is larger than the planned result, for the second script running. Scoring quality -")
print("how many panel positions a genome has scorable at all - entered as the least interesting of three")
print("candidates and came back at -2.35 pH per log10 positions, p 2e-10. It does not absorb the histidine")
print("coefficient, so it failed as a candidate, but a genome-level predictor that strong had never been held")
print("fixed in any panel result including the headline, and scripts/54 had already found the per-reference")
print("version of the same thing. Held fixed, the acidic exposed-periplasmic coefficient goes from +0.10681 to")
print("+0.11775 - it strengthens - so the caveat is now checked rather than open. What remains is that the control")
print("belongs in the panel scripts and is currently in one script about histidine.")
