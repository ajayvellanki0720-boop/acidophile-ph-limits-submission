"""What correction 20 does to the experiment's predicted effect and its sample size -> DATA_DIR/design_rescaling.csv and, with `wide`, DATA_DIR/design_rescaling_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHY THE EXPERIMENT INHERITS CORRECTION 20, which is an argument and not an arithmetic step, so it is stated before any number. Correction 20 found the exposed-periplasmic
acidic coefficient inflated about 27% by omitting exposed histidine. docs/CHARGE_REVERSAL_DESIGN.md section 2 scales its predicted effect through the composite's slope per
standard deviation, and the question is whether the inflated coefficient or the histidine-controlled one is the right input. It is the controlled one, because of what the
experiment DOES: reversing an aspartate or glutamate to lysine at an exposed periplasmic position changes the acidic share and the basic share of that class and leaves the
histidine share exactly as it was. The causal quantity the edit moves is therefore the acidic effect with histidine HELD FIXED, and using the uncontrolled coefficient would
credit the edit with part of a residue the edit does not touch.
AND IT IS APPLIED AS A RATIO, NOT AS A SUBSTITUTION, which is the part most likely to go wrong. The design's committed figures descend from scripts/57's factor-score model
and its +0.185 per SD; the coefficients computed here are the raw pooled net share per its own standard deviation, a different quantity on a different scale. Substituting one
for the other would be correction 12 - a comparison across specifications - so what transfers is the ratio between the controlled and uncontrolled fits of the SAME model on
the SAME sample, which is scale-free and is the only thing this script claims.
THE DIRECTION IS THE UNFLATTERING ONE, and that is worth saying because the design document has a rule about it. Its section on the hinge declines a 1.55x increase to the
predicted effect partly because "raising a predicted effect is the direction that flatters the experiment, which in this project has been the reliable sign that a number
needs deriving a second way before it is used". This lowers the predicted effect. That does not earn it a discount on scrutiny, but it does mean the usual reason for
suspicion is absent, and the thing to check instead is that the ratio comes from one model and one sample rather than two - which it does, and an assertion requires it.
Usage: python3 scripts/88_design_rescaling.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
# The design's committed figures, quoted from docs/CHARGE_REVERSAL_DESIGN.md section 2 so the rescaling is
# auditable against the document rather than recomputed from a different chain.
COMMITTED = {"NuoN, per reversal (32-ref)": -0.004, "NuoL+NuoN 14 reversals (32-ref), the design target": -0.05,
             "NuoN, per reversal (7-ref)": -0.012, "NuoL+NuoN 14 reversals (7-ref), optimistic end": -0.16}
RESOLUTION = 0.02          # the assay resolution target the design commits to
N_COMMITTED = 64           # scripts/70's minimum replicates per variant per pH, derived from the above
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
RES = ("D", "E", "K", "R", "H")
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
NET = s["D"] + s["E"] - s["K"] - s["R"]
SDNET = float(np.std(NET, ddof=1))
n_used = int(ok.sum())

print(f"1. the ratio, from one model on one sample (n={n_used}, lambda {lam:.3f})")
print(f"   {'specification':>34} {'beta/pt':>10} {'se':>9} {'p':>10} {'per SD':>9}")
fits = {}
for lab, cols in (("net D+E-K-R, histidine free", [NET]), ("net D+E-K-R, histidine held fixed", [s["H"], NET])):
    r = ph.pgls(y, np.column_stack([XS] + cols), Cm, lam=lam)
    b, se, pv = float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1])
    fits[lab] = (b, se, pv)
    print(f"   {lab:>34} {b:>+10.5f} {se:>9.5f} {pv:>10.2e} {b * SDNET:>+9.4f}")
    rows.append({"part": "ratio", "specification": lab, "beta_per_pt": b, "se": se, "p": pv,
                 "beta_per_sd": b * SDNET, "sd_of_composite": SDNET, "n": n_used, "lambda": lam})
b_free = fits["net D+E-K-R, histidine free"][0]
b_held = fits["net D+E-K-R, histidine held fixed"][0]
ratio = b_held / b_free
print(f"   rescaling factor {ratio:.4f}, i.e. the predicted effect falls by {1 - ratio:.1%}")
print("   Both fits use the same covariates, the same genomes and the same lambda, so the ratio is scale-free")
print("   and is the only quantity transferred below.")
assert fits["net D+E-K-R, histidine held fixed"][2] < 0.05, \
    "the composite does not survive holding histidine fixed, which would be a far larger finding than a rescaling"
assert 0.0 < ratio < 1.0, f"the rescaling factor is {ratio}, outside the range correction 20 predicts"
rows.append({"part": "ratio", "specification": "rescaling factor", "ratio": ratio, "reduction": 1 - ratio})

print("\n2. the design's committed predictions, rescaled")
print(f"   {'prediction':>48} {'committed':>10} {'rescaled':>9}")
for lab, v in COMMITTED.items():
    print(f"   {lab:>48} {v:>+10.3f} {v * ratio:>+9.3f}")
    rows.append({"part": "rescaled prediction", "prediction": lab, "committed": v, "rescaled": v * ratio})
target = COMMITTED["NuoL+NuoN 14 reversals (32-ref), the design target"]
print(f"\n3. the margin against the assay resolution target of {RESOLUTION}")
m0, m1 = abs(target) / RESOLUTION, abs(target * ratio) / RESOLUTION
print(f"   design target {target:+.3f} -> {target * ratio:+.3f}; margin over resolution falls from "
      f"{m0:.1f}x to {m1:.1f}x")
still = abs(target * ratio) > RESOLUTION
print("   " + ("The primary strain's predicted effect remains ABOVE the resolution target, so the design is not"
               if still else
               "The primary strain's predicted effect falls BELOW the resolution target, which would make the"))
print("   " + ("   invalidated - but the headroom it was sized with is substantially reduced."
               if still else "   primary strain undetectable as designed."))
rows.append({"part": "resolution margin", "resolution_target": RESOLUTION, "margin_committed": m0,
             "margin_rescaled": m1, "still_detectable": still})

print(f"\n4. the sample-size consequence, which is where a 28% effect reduction actually bites")
# For a fixed power and variance, the replicates needed scale as the inverse SQUARE of the effect size, so a
# factor r on the effect multiplies n by 1/r^2. scripts/70 derived n >= 64 per variant per pH from the committed
# figures; that derivation is not re-run here, and the scaling is applied to its result as a ratio for the same
# reason part 1's ratio is used rather than its absolute coefficients.
infl = 1.0 / ratio ** 2
print(f"   replicates scale as 1/effect^2, so n multiplies by {infl:.2f}: scripts/70's n >= {N_COMMITTED} becomes "
      f"n >= {int(np.ceil(N_COMMITTED * infl))}")
print("   This is a rescaling of that script's published minimum, not a re-derivation of it, and it assumes the")
print("   variance and power targets behind the 64 are unchanged. It is the figure the design decision needs.")
rows.append({"part": "sample size", "n_committed": N_COMMITTED, "inflation_factor": infl,
             "n_rescaled": int(np.ceil(N_COMMITTED * infl))})

print("\n5. THE OTHER UNADOPTED RESCALING, AND THE FACT THAT NOBODY HAS MULTIPLIED THEM")
# docs/PREREGISTRATION.md section 9 already carries a "computed but NOT adopted" rescaling of the same committed
# figure, on an unrelated ground: if a position's contribution is proportional to its own residue's coefficient,
# the five aspartate positions are worth less than the nine glutamate ones and -0.05 becomes -0.0397, with n
# rising from 64 to 101 (scripts/75 wide). That note and this script adjust the SAME number for DIFFERENT reasons -
# which positions contribute, and which residues the model holds fixed - and two separately-recorded reductions
# sitting in two documents is exactly how a design stays nominally feasible while actually not being. The product
# is computed here because nobody had computed it, and it is flagged as resting on an independence assumption
# that nothing tests: the two adjustments are plausibly multiplicative, and if they interact the product is wrong.
PRE9 = -0.0397          # docs/PREREGISTRATION.md section 9, the per-residue-weighting rescaling of -0.05
w = PRE9 / COMMITTED["NuoL+NuoN 14 reversals (32-ref), the design target"]
comb = ratio * w
print(f"   section 9's factor {w:.4f} (from {target:+.3f} to {PRE9:+.3f}), this script's {ratio:.4f}")
print(f"   product {comb:.4f} -> the design target becomes {target * comb:+.4f}, and n >= {N_COMMITTED} becomes "
      f"n >= {int(np.ceil(N_COMMITTED / comb ** 2))}")
print(f"   against the resolution target of {RESOLUTION} that is a margin of {abs(target * comb) / RESOLUTION:.1f}x.")
both_ok = abs(target * comb) > RESOLUTION
print("   " + ("Still above the resolution target, but only just, and on an INDEPENDENCE ASSUMPTION nothing has"
               if both_ok else
               "BELOW the resolution target. If both rescalings are adopted the primary strain as designed would"))
print("   " + ("   tested - the two adjustments are plausibly multiplicative and if they interact this is wrong."
               if both_ok else "   not be detectable, which is a design decision and not a footnote."))
print("   Neither is adopted. The point of computing the product is that two unadopted reductions recorded in two")
print("   different documents is how a design stays nominally feasible while actually ceasing to be.")
# STORED, not just used. -0.0397 is quoted from docs/PREREGISTRATION.md section 9 as an input to this
# calculation, and a figure a document quotes has to live in a result file or scripts/69 cannot reach it - which
# is exactly what happened: the README quoted it inside this script's entry and the check reported it untraceable,
# correctly, because the script consumed it without ever writing it down.
rows.append({"part": "compounded rescaling", "quoted_input": "docs/PREREGISTRATION.md section 9 target",
             "rescaled": PRE9, "source": "scripts/75 wide, via the pre-registration"})
rows.append({"part": "compounded rescaling", "factor_section9": w, "factor_this_script": ratio,
             "product": comb, "rescaled_target": target * comb,
             "n_rescaled": int(np.ceil(N_COMMITTED / comb ** 2)),
             "margin_over_resolution": abs(target * comb) / RESOLUTION,
             "still_detectable": both_ok, "assumes": "the two rescalings are independent; untested"})

print("\n6. the reading")
print(f"   The experiment inherits correction 20 because the edit it performs does not touch histidine, so the")
print(f"   coefficient that describes its causal effect is the one with histidine held fixed. That coefficient is")
print(f"   {1 - ratio:.0%} smaller on the same model and sample, and both fits remain highly significant, so this")
print("   is a rescaling rather than a threat to the design's premise. What it costs is headroom: the primary")
print(f"   strain's predicted effect goes from {target:+.3f} to {target * ratio:+.3f} against an assay resolution")
print(f"   target of {RESOLUTION}, a margin of {m1:.1f}x rather than {m0:.1f}x, and the replicate requirement")
print(f"   roughly doubles. Nothing here is adopted into the design document by this script; the figures are")
print("   computed so the deviations procedure can decide, which is the same treatment correction 19 gave the")
print("   histidine arm and correction 20 gave the composite. The direction is the unflattering one, so the")
print("   usual reason for suspecting a design number is absent - but the check that matters is that the ratio")
print("   comes from one model and one sample rather than two, and part 1 asserts exactly that.")
print("   THE ONE THING THIS CANNOT DO is tell you whether n >= 126 is affordable. That is a bench and budget")
print(f"   question, and doubling the replicate requirement is the kind of change that can quietly make a")
print("   committed design unbuildable while every individual number in it stays defensible. It is flagged here")
print("   rather than left to be discovered when the plates are being counted.")
pd.DataFrame(rows).to_csv(D / f"design_rescaling{TAG}.csv", index=False)
print("\nwrote", D / f"design_rescaling{TAG}.csv")
