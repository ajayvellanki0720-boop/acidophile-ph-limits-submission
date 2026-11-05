"""The pKa test, run POOLED after the per-reference version refused itself -> DATA_DIR/histidine_pooled.csv and, with the `wide` argument, DATA_DIR/histidine_pooled_wide.csv
BOTH SPELLINGS ARE DECLARED because scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate from the default run's. Declaring only the default is what made the README's three correction-18 coefficients untraceable the moment they were quoted: the figures were right, the file they
live in is the _wide one, and the pool built from the base spelling was empty. scripts/77 check C4 had this flagged as latent debt beforehand - latent debt is exactly
what it was, and quoting from the file is what turned it into a miss.
THIS CORRECTS A DIAGNOSIS MADE EARLIER THE SAME DAY, not a number. scripts/76 declared a floor before the data existed - the median reference must carry at least 2 histidines
at exposed periplasmic positions - and refused the test when the observed median came back 0.00. The floor was right: nothing can be estimated per reference off counts of
0 or 1. The DIAGNOSIS attached to it was wrong. It was recorded as "histidine is essentially absent from these surfaces", and as a statement about composition that is false:
histidine is 1.69% of exposed periplasmic positions, HIGHER than its 1.22% on the surfaces of the soluble control proteins. The median is zero because each reference
contributes only about 27 exposed periplasmic positions, so about 0.46 histidines are expected per reference per genome. The constraint was class SIZE, never absence.
AND THE REMEDY IS THE ONE THIS PROJECT ALREADY USES EVERYWHERE ELSE. The headline does not read a residue share off one reference; it pools 32. Pooling the histidine count
and the position count across references before taking the share gives roughly 860 positions per genome, which is estimable. So the question scripts/76 left open is
answerable, and it was closed for the wrong reason - a refusal is only as good as the account of what is blocking it.
THE PREDICTION WAS FIXED IN ADVANCE AND IS GENUINELY TWO-SIDED, which is why the test discriminates instead of confirming. It is recorded in scripts/36a at the tally and in
docs/SYNTHESIS.md section 8. Side-chain pKa: D ~3.65, E ~4.25, H ~6.0. Glutamate carries about 2.4x what aspartate does (scripts/75), which is what a CHARGE-BURDEN account
predicts - the carboxylate whose charge state turns over first as the envelope acidifies is the one that matters - and on that account histidine, turning over earlier still,
should behave like a stronger glutamate: SAME SIGN, larger. A BUFFERING account predicts the OPPOSITE sign, since histidine is the main protein buffer in this pH range and
more of it should help rather than hurt. This project's own withdrawn histidine lead ran the buffering way proteome-wide. Either outcome is informative; neither is the
outcome the project would prefer.
SCALE. Every coefficient is per percentage point, from ph.pgls directly rather than ph.pgls_scan, so correction 15's mislabelling cannot arise. D, E and H are entered in the
SAME model so each is that residue's effect holding the others fixed, which is correction 14's lesson: aspartate looked like a coin flip only because glutamate was omitted.
Usage: python3 scripts/81_histidine_pooled.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
# The histidine column exists only in the rebuilt panel. The committed panel predates it, so this reads the
# DETERMINISTIC build (correction 17) - which is also the only panel that is a function of its inputs.
PANEL = "structure_panel_wide_det.csv" if WHICH == "wide" else "structure_panel.csv"
TAG = "_wide" if WHICH == "wide" else ""
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
SOL = [str(x) for x in _r[_r.accepted & ~_r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
if not any(c.endswith("|H") for c in pan.columns):
    sys.exit(f"{PANEL} has no histidine columns; rebuild the panel with scripts/36a first")

print("1. the compositional fact that retires \"histidine is absent from these surfaces\"")
print(f"   {'class':>10} {'refs':>5} {'positions':>12} {'histidines':>11} {'share %':>9} {'per ref/genome':>15}")
for cls in ("peri_exp", "cyto_exp", "peri_bur", "cyto_bur", "TM", "surface", "buried"):
    H = N = 0.0
    k = 0
    for ref in MEM + SOL:
        hc, nc = f"{ref}|{cls}|H", f"{ref}|{cls}|n"
        if hc in pan and nc in pan:
            H += pan[hc].sum(); N += pan[nc].sum(); k += 1
    if not k:
        continue
    print(f"   {cls:>10} {k:>5} {N:>12.0f} {H:>11.0f} {100 * H / N:>9.3f} {H / N * (N / k / len(pan)):>15.2f}")
    rows.append({"part": "composition", "class": cls, "refs": k, "positions": N, "histidines": H,
                 "share_pct": 100 * H / N})
pe = [r for r in rows if r["class"] == "peri_exp"][0]
sf = [r for r in rows if r["class"] == "surface"]
print(f"   Exposed periplasmic is {pe['share_pct']:.2f}% histidine"
      + (f", against {sf[0]['share_pct']:.2f}% on the soluble controls' surfaces." if sf else "."))
print("   So the per-reference median of zero was CLASS SIZE, not absence: about "
      f"{pe['positions'] / pe['refs'] / len(pan):.0f} positions per reference per genome at that share is "
      f"{pe['share_pct'] / 100 * pe['positions'] / pe['refs'] / len(pan):.2f} histidines.")

print("\n2. the pooled test, D/E/H entered together (per percentage point, composition and source controlled)")
df = (pd.read_csv(D / "phenotypes_all.csv")
      # frac_His is merged because part 4 needs the genome's OWN histidine share: `gw` is genome-wide ACIDIC
      # excess and controls nothing about histidine, and this project has already withdrawn a proteome-wide
      # histidine lead, so that is the named confound rather than a generic worry.
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp",
                                                     "frac_His"]], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]

# POOLED SHARES. Sum the counts and the positions across references, then divide - not the mean of per-reference
# shares, which weights a reference with 9 scorable positions like one with 80 and is what made the per-reference
# version unusable. A genome needs a minimum pooled position count to be scored at all, for the same reason.
MINPOS = 200
tot = {k: np.zeros(len(B)) for k in ("n", "D", "E", "H", "K")}
used = 0
for ref in MEM:
    nc = f"{ref}|peri_exp|n"
    if nc not in B or f"{ref}|peri_exp|H" not in B:
        continue
    n = B[nc].fillna(0).to_numpy(float)
    tot["n"] += n
    for k in ("D", "E", "H", "K"):
        tot[k] += B[f"{ref}|peri_exp|{k}"].fillna(0).to_numpy(float)
    used += 1
ok = tot["n"] >= MINPOS
print(f"   {used} references pooled; {int(ok.sum())} of {len(B)} genomes reach {MINPOS} pooled positions")
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
X0 = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
sh = {k: 100 * tot[k][ok] / tot["n"][ok] for k in ("D", "E", "H", "K")}
y = g.ph_min.to_numpy(float)
lam = ph.pgls(y, X0, Cm)["lambda"]
r = ph.pgls(y, np.column_stack([X0, sh["D"], sh["E"], sh["H"], sh["K"]]), Cm, lam=lam)
print(f"   {'residue':>9} {'pKa':>5} {'share %':>9} {'beta/pt':>10} {'se':>9} {'p':>10}")
PKA = {"D": 3.65, "E": 4.25, "H": 6.0, "K": 10.5}
est = {}
for i, k in enumerate(("D", "E", "H", "K")):
    b, se, p = float(r["beta"][-4 + i]), float(r["se"][-4 + i]), float(r["p"][-4 + i])
    est[k] = (b, se, p)
    print(f"   {k:>9} {PKA[k]:>5.2f} {sh[k].mean():>9.3f} {b:>+10.5f} {se:>9.5f} {p:>10.2e}")
    rows.append({"part": "pooled joint fit", "residue": k, "pka": PKA[k], "mean_share_pct": float(sh[k].mean()),
                 "beta_per_pt": b, "se": se, "p": p, "n": len(g), "lambda": lam})

print("\n3. which account the sign picks out, against the prediction fixed in advance")
bH, seH, pH = est["H"]
bE = est["E"][0]
same = (bH > 0) == (bE > 0)
if pH >= 0.05:
    verdict = ("INDETERMINATE: histidine's coefficient is not distinguishable from zero, so neither account is "
               "supported. The test ran, and it does not decide.")
elif same and abs(bH) > abs(bE):
    verdict = ("CHARGE BURDEN, in the strong form: same sign as glutamate and larger, which is what a residue "
               "titrating earlier should do if the mechanism is charge.")
elif same:
    verdict = ("CHARGE BURDEN, in the weak form: same sign as glutamate but not larger, so the sign follows "
               "charge while the pKa ORDERING does not.")
else:
    verdict = ("BUFFERING: opposite sign to glutamate, so more histidine goes with a lower pH limit rather than "
               "a higher one - the direction this project's withdrawn proteome-wide lead ran.")
print(f"   glutamate {bE:+.5f}/pt, histidine {bH:+.5f}/pt (p {pH:.2e})")
print(f"   {verdict}")
rows.append({"part": "verdict", "beta_E_per_pt": bE, "beta_H_per_pt": bH, "p_H": pH,
             "same_sign_as_glutamate": same, "verdict": verdict})

print("\n4. THE CONTROL THAT DECIDES IT: genome-wide histidine")
# Part 2's covariate `gw` is genome-wide ACIDIC excess. It does not control genome-wide HISTIDINE at all, so the
# periplasmic histidine share partly measures how much histidine the genome has anywhere. This project has already
# withdrawn a proteome-wide histidine lead, so that is the specific confound to rule out rather than a generic
# worry - and if the surface share survives controlling the genome's own share, the effect is about the surface.
gwh = g.frac_His.to_numpy(float) * 100
r2 = ph.pgls(y, np.column_stack([X0, gwh, sh["D"], sh["E"], sh["H"]]), Cm, lam=lam)
b2, se2, p2 = float(r2["beta"][-1]), float(r2["se"][-1]), float(r2["p"][-1])
bg, pg = float(r2["beta"][-4]), float(r2["p"][-4])
rr = float(np.corrcoef(gwh, sh["H"])[0, 1])
print(f"   correlation between the genome's histidine share and the surface share: r {rr:+.3f}")
print(f"   genome-wide histidine itself:   {bg:+.5f}/pt (p {pg:.2e})")
print(f"   surface histidine, controlled:  {b2:+.5f}/pt (se {se2:.5f}, p {p2:.2e})"
      f"   [uncontrolled {bH:+.5f}]")
surv = p2 < 0.05 and (b2 > 0) == (bH > 0)
print(f"   {'SURVIVES' if surv else 'DOES NOT SURVIVE'}: "
      + (f"{abs(b2 / bH):.0%} of the uncontrolled coefficient" if surv else
         "the surface share carries no information once the genome's own histidine content is controlled"))
rows.append({"part": "genome-wide histidine control", "beta_H_controlled_per_pt": b2, "se": se2, "p": p2,
             "beta_genomewide_His_per_pt": bg, "p_genomewide": pg, "r_surface_vs_genomewide": rr,
             "survives": surv, "share_of_uncontrolled": abs(b2 / bH)})

print("\n5. the side-of-membrane control, which asks whether this is specific at all")
# A genome-wide histidine trend would show up in EVERY class. The exposed CYTOPLASMIC class is the project's
# standing control: same proteins, same alignments, same scoring, opposite side of the membrane.
tc = {k: np.zeros(len(B)) for k in ("n", "H")}
for ref in MEM:
    if f"{ref}|cyto_exp|H" not in B:
        continue
    tc["n"] += B[f"{ref}|cyto_exp|n"].fillna(0).to_numpy(float)
    tc["H"] += B[f"{ref}|cyto_exp|H"].fillna(0).to_numpy(float)
okc = (tc["n"] >= MINPOS) & ok
gc_ = B[okc].reset_index(drop=True)
Cc = ph.vcv_fast(tree, list(gc_.ncbi_assembly))
Xc = np.column_stack([np.ones(len(gc_))] + [gc_[k] for k in ["gw"] + SRCC])
hc = 100 * tc["H"][okc] / tc["n"][okc]
lamc = ph.pgls(gc_.ph_min.to_numpy(float), Xc, Cc)["lambda"]
rc = ph.pgls(gc_.ph_min.to_numpy(float), np.column_stack([Xc, hc]), Cc, lam=lamc)
bc, sec, pc = float(rc["beta"][-1]), float(rc["se"][-1]), float(rc["p"][-1])
print(f"   exposed cytoplasmic histidine: {bc:+.5f}/pt (se {sec:.5f}, p {pc:.2e}), n={len(gc_)}")
print(f"   periplasmic against cytoplasmic: {bH:+.5f} against {bc:+.5f}")
rows.append({"part": "cytoplasmic control", "beta_per_pt": bc, "se": sec, "p": pc, "n": len(gc_)})

# SPECIFICITY, judged by the project's own standard. The acidic finding's exposed-cytoplasmic control is null
# (+0.019, p 0.38). If histidine's control is significant and same-signed, the association is not specific to the
# outward-facing surface, whatever its p-value on that surface - and a non-specific surface effect is what a
# genome-wide compositional trend looks like after it has been sliced by position class.
spec = not (pc < 0.05 and (bc > 0) == (bH > 0))
print(f"   specificity by the standard applied to the acidic result: "
      + ("PASSES - the cytoplasmic control is null" if spec else
         f"FAILS - the control is significant (p {pc:.1e}) and the same sign, at {abs(bc / bH):.0%} of the "
         f"periplasmic coefficient"))
rows.append({"part": "specificity", "passes": spec, "control_share_of_periplasmic": abs(bc / bH)})

print("\n6. reconciling the D and E coefficients with scripts/75, which look inconsistent and are not")
# scripts/75 pools PER-REFERENCE fits and reports D +0.00500 and E +0.01185 per point on this panel; the fit above
# uses one POOLED share and gives figures several times larger. That is not a contradiction, it is correction 4's
# finding restated: a per-reference share is a noisy measure of the same quantity, and classical attenuation
# shrinks its coefficient towards zero. The invariant to check is the RATIO, which the scale cancels out of.
R75 = {"D": 0.00500, "E": 0.01185}
r75 = R75["E"] / R75["D"]
rhere = est["E"][0] / est["D"][0]
print(f"   per-reference pooling (scripts/75): E/D = {r75:.2f}x     pooled share (here): E/D = {rhere:.2f}x")
print(f"   the two agree to {abs(rhere - r75) / r75:.0%}, so the glutamate-over-aspartate ranking replicates")
print("   under a different estimator; only the scale differs, as measurement error predicts.")
rows.append({"part": "reconciliation with scripts/75", "ratio_E_over_D_per_reference_pooled": r75,
             "ratio_E_over_D_pooled_share": rhere, "relative_difference": abs(rhere - r75) / r75})

print("\n7. the design implication, COMPUTED AND DELIBERATELY NOT ADOPTED")
# The experiment's primary arm substitutes acidic -> LYSINE, and correction 14 already cost that choice half its
# rationale by showing lysine does not replicate on this panel. Histidine is the obvious alternative basic-ish
# substitution, and on THIS estimator - the same model, the same sample, so correction 12's across-sample error
# cannot arise - it carries about ten times lysine's coefficient and is significant where lysine's is not.
bK, seK, pK = est["K"]
print(f"   lysine   {bK:+.5f}/pt (p {pK:.2f}){'  not significant' if pK >= 0.05 else ''}")
print(f"   histidine {bH:+.5f}/pt (p {pH:.1e})   -> {abs(bH / bK):.1f}x lysine's magnitude")
print("   Both are negative, i.e. both go with a LOWER pH limit, which is the direction the reversal wants, so")
print("   a D/E->H arm would be predicted to shift further than the committed D/E->K arm.")
print("   IT IS NOT ADOPTED, for a reason that is not caution for its own sake: histidine FAILS the specificity")
print("   control that the acidic result passes, so part of its coefficient is not attributable to the outward-")
print("   facing surface the experiment edits, and substituting on a coefficient that is partly non-specific")
print("   would be building the arm on the weakest number in this script rather than the strongest. Recorded")
print("   here so it is a decision available to the deviations procedure rather than something to rediscover at")
print("   the bench - the same treatment as the n 64->101 note in docs/PREREGISTRATION.md section 9.")
rows.append({"part": "design implication, not adopted", "beta_K_per_pt": bK, "p_K": pK,
             "beta_H_per_pt": bH, "p_H": pH, "H_over_K_magnitude": abs(bH / bK),
             "adopted": False, "why_not": "histidine fails the exposed-cytoplasmic specificity control"})

pd.DataFrame(rows).to_csv(D / f"histidine_pooled{TAG}.csv", index=False)
print("\nwrote", D / f"histidine_pooled{TAG}.csv")
print("\nReading, and it is a mixed verdict that must not be rounded to either edge.")
print("scripts/76 was right to refuse the per-reference test and wrong about why. \"Absent\" closes a question;")
print("\"too few positions per reference\" names a fix the project already uses for every other residue, and")
print("pooling is that fix. Pooled, the test runs and the pre-registered pKa ordering FAILS: charge burden")
print("predicted histidine would behave like a stronger glutamate - same sign, larger - and it is the opposite")
print("sign at twice the magnitude. The buffering direction is what the data show. Two things keep that from")
print("being a mechanism claim. It is NOT the project's withdrawn proteome-wide histidine lead, which is the")
print("one confound worth naming here: the genome's own histidine share predicts nothing (p 0.96) and the")
print("surface share survives controlling it at 100% of its uncontrolled size. But it FAILS the specificity")
print("standard the acidic result passes - 55% of the coefficient appears on the cytoplasmic side too, where")
print("the acidic control is null. A charge-burden mechanism should be strictly periplasmic, since that is")
print("where the proton gradient is; buffering need not be. So the sign and the partial two-sidedness point the")
print("same way, which is suggestive and is not proof, and the right tier for this is a real association whose")
print("mechanism is unsettled - not a confirmed mechanism, and not a null. What it does settle is that the")
print("residue ranking does NOT follow side-chain pKa monotonically, which was the question asked.")
print("The glutamate-over-aspartate ranking meanwhile replicates under this estimator, at a ratio within 10%")
print("of scripts/75's, which is the useful by-product of having had to build a second route to the same thing.")
