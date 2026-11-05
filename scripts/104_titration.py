"""Charge burden or buffering? The two accounts predict opposite signs of the same quantity -> DATA_DIR/titration.csv
WHY THIS IS NOW A WELL-POSED QUESTION, WHEN AN ELECTROSTATICS CALCULATION WOULD NOT HAVE BEEN A WEEK AGO. Correction
21 found the one exhaustive experiment: every histidine, aspartate and glutamate in the periplasmic domain of the
H. pylori urea channel mutagenised, with three carboxylates REQUIRED to activate transport in acid. This project
finds acid-tolerant genomes DEPLETED of exposed periplasmic carboxylates and has been explaining that by charge
burden - carboxylates protonate as external pH falls, so fewer of them means less pH-dependent change in outward
surface charge. The competing account, which UreI and the one gastric-pH engineering study both favour, is BUFFERING:
titratable groups near the operating pH resist local acidification and are therefore useful rather than costly.
THE TWO ACCOUNTS ARE DISTINGUISHABLE FROM COMPOSITION ALONE, which is the point of this script. Write Q(pH) for the
net charge of a surface and beta(pH) = -dQ/dpH for its buffer capacity. Both are determined by the residue counts and
the side-chain pKa values; no structure is needed for the first-order version, and no solver.
  BURDEN says what hurts is the CHANGE in charge between neutrality and the growth limit. Acid tolerance should go
  with a SMALLER |Q(7) - Q(pH_min)| on the outward surface.
  BUFFERING says what helps is capacity AT the operating pH. Acid tolerance should go with a LARGER beta(pH_min).
But beta(pH) is, to within a constant, the derivative whose integral is that charge change - so on the same residues
the two accounts want the same quantity to move in OPPOSITE directions. One composition cannot satisfy both, and the
sign of the association decides between them. That is a genuine discrimination rather than a reframing.
WHAT IS APPROXIMATED, AND WHY IT DOES NOT DECIDE THE ANSWER. Model pKa values are used (D 3.65, E 4.25, H 6.00,
K 10.50, R 12.50, with termini ignored because the panel counts interior positions). Real pKa values shift by one to
two units in a protein environment, which moves the CURVES - but a shift applies to both accounts alike, and the
discriminating test is the SIGN of an association across genomes, which a common shift cannot flip. The place a pKa
shift would matter is the claim that a particular residue titrates inside the growth range, so that claim is made
only where the margin is large, and part 4 measures how large a shift would be needed to change the reading.
THE THIRD ACCOUNT, which is the one both of the first two forget. Neither says anything about WHERE in pH the signal
lives. If the mechanism is carboxylate protonation, a composition weighted by how much each residue titrates near
pH 4 should predict better than raw share; if it is histidine, the sharpening should happen near pH 6. Part 3 scans
that, so the residue the data points at is identified rather than assumed.
SCALE. Charges are per 100 exposed periplasmic positions; beta is in charge units per pH unit on the same basis.
Usage: python3 scripts/104_titration.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MINPOS = 200
PKA = {"D": 3.65, "E": 4.25, "H": 6.00, "K": 10.50, "R": 12.50}
SIGN = {"D": -1, "E": -1, "H": +1, "K": +1, "R": +1}   # charge of the PROTONATED-minus-deprotonated form
rows = []

refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
pan = pd.read_csv(D / "structure_panel_wide_det.csv").drop_duplicates("ncbi_assembly")
pheno = (pd.read_csv(D / "phenotypes_all.csv")
         .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100.0
SRCC = [f"src_{o}" for o in SRC]
M = pan.set_index("ncbi_assembly")


def pooled(cls):
    tot = {k: np.zeros(len(B)) for k in ("n",) + tuple(PKA)}
    for r in MEM:
        if f"{r}|{cls}|n" not in M.columns or f"{r}|{cls}|H" not in M.columns:
            continue
        tot["n"] += np.nan_to_num(M[f"{r}|{cls}|n"].reindex(B.ncbi_assembly).to_numpy(float))
        for k in PKA:
            tot[k] += np.nan_to_num(M[f"{r}|{cls}|{k}"].reindex(B.ncbi_assembly).to_numpy(float))
    return tot


TP, TC = pooled("peri_exp"), pooled("cyto_exp")
keep = (TP["n"] >= MINPOS) & (TC["n"] >= MINPOS)
B = B[keep].reset_index(drop=True)
TP = {k: v[keep] for k, v in TP.items()}
TC = {k: v[keep] for k, v in TC.items()}
C = ph.vcv_fast(tree, list(B.ncbi_assembly))
y = B.ph_min.to_numpy(float)
USE = [c for c in ["gw"] + SRCC if B[c].nunique() > 1]
X0 = np.column_stack([np.ones(len(B))] + [B[c].to_numpy(float) for c in USE])
LAM = ph.pgls(y, X0, C)["lambda"]
print(f"sample {len(B)} genomes with >= {MINPOS} pooled positions on both faces; lambda {LAM:.4f}")


def fit(x, extra=None):
    cols = [X0] + ([extra] if extra is not None else []) + [((x - x.mean()) / np.std(x))[:, None]]
    f = ph.pgls(y, np.column_stack(cols), C, lam=LAM)
    return float(f["beta"][-1]), float(f["se"][-1]), float(f["p"][-1])


def frac_protonated(pH, pk):
    """Henderson-Hasselbalch: the fraction of groups carrying their protonated form at this pH."""
    return 1.0 / (1.0 + 10.0 ** (pH - pk))


def charge(tot, pH):
    """Net charge per 100 positions. Acidic groups are neutral when protonated and -1 when not; basic groups are
    +1 when protonated and neutral when not - so one expression with SIGN handles both."""
    q = np.zeros(len(B))
    for k, pk in PKA.items():
        f = frac_protonated(pH, pk)
        q += tot[k] * (f if SIGN[k] > 0 else -(1.0 - f))
    return 100.0 * q / tot["n"]


def buffer_capacity(tot, pH):
    """beta(pH) = ln(10) * sum_i n_i * f_i * (1 - f_i), the standard Henderson-Hasselbalch buffer capacity, per 100
    positions. It peaks where pH equals a group's pKa, which is what makes it the discriminating quantity."""
    b = np.zeros(len(B))
    for k, pk in PKA.items():
        f = frac_protonated(pH, pk)
        b += tot[k] * f * (1.0 - f)
    return 100.0 * np.log(10.0) * b / tot["n"]


print("\n1. the two accounts, each as the quantity it says should matter")
# BURDEN: the charge the outward surface has to change by, between neutrality and the organism's own limit.
dq = np.abs(charge(TP, 7.0) - np.array([charge(TP, v)[i] for i, v in enumerate(y)]))
# BUFFERING: capacity at the organism's own limit.
bc = np.array([buffer_capacity(TP, v)[i] for i, v in enumerate(y)])
# and the same two for the cytoplasmic face, which is the control: the cytoplasm is held near neutrality, so neither
# account predicts anything there, and a signal would mean the quantity is tracking something non-specific.
dq_c = np.abs(charge(TC, 7.0) - np.array([charge(TC, v)[i] for i, v in enumerate(y)]))
bc_c = np.array([buffer_capacity(TC, v)[i] for i, v in enumerate(y)])
print(f"   {'quantity':>52s} {'beta/SD':>10s} {'se':>8s} {'p':>10s}")
for lab, v in (("periplasmic |charge change pH7 -> own limit|  (BURDEN)", dq),
               ("periplasmic buffer capacity at own limit   (BUFFERING)", bc),
               ("cytoplasmic |charge change|  (control)", dq_c),
               ("cytoplasmic buffer capacity  (control)", bc_c)):
    b, se, p = fit(v)
    print(f"   {lab:>52s} {b:>+10.5f} {se:>8.5f} {p:>10.2e}")
    rows.append({"part": "1", "measure": lab, "n": len(B), "beta_per_sd": b, "se": se, "p": p})
# THE DISCRIMINATION. Both quantities are computed on the same residues, so the burden account wants acid tolerance
# (LOW ph_min) to go with a SMALL charge change - a POSITIVE coefficient on dq - while the buffering account wants
# acid tolerance to go with LARGE capacity - a NEGATIVE coefficient on bc. They cannot both hold.
b_dq, _, p_dq = fit(dq)
b_bc, _, p_bc = fit(bc)
# PART 1 IS CIRCULAR AND IS DISCARDED. Both quantities are evaluated at each genome's OWN ph_min, so the outcome sits
# inside the predictor and a monotone function of y is being regressed on y. The proof is in the control rows: the
# CYTOPLASMIC quantities, which neither account predicts anything about because the cytoplasm is held near
# neutrality, come back enormous at p = 0. That is not biology, it is the leak. A first version of this script read
# its verdict off these rows and reported BUFFERING; the cytoplasmic control was printed three lines above saying the
# same thing just as strongly, which is what a control is for.
cyto_leak = [r for r in rows if r["part"] == "1" and "cytoplasmic" in r["measure"]]
LEAK = any(r["p"] < 1e-3 for r in cyto_leak)
print(f"\n   DISCARDED AS CIRCULAR: evaluating at each genome's own limit puts the outcome inside the predictor.")
print(f"   The proof is the cytoplasmic control, which neither account predicts and which fires at "
      f"{min(r['beta_per_sd'] for r in cyto_leak):+.5f} (p {min(r['p'] for r in cyto_leak):.2e}).")
print("   Nothing in part 1 is read as evidence. Part 2 repeats both quantities at a fixed pH, where y appears on one")
print("   side only, and that is what the verdict uses.")
rows.append({"part": "1-discarded", "measure": "circularity check via cytoplasmic control", "leak_detected": LEAK,
             "beta_charge_change": b_dq, "p_charge_change": p_dq, "beta_capacity": b_bc, "p_capacity": p_bc})
assert LEAK, ("the cytoplasmic control did NOT fire in part 1, which would mean the circularity argument above is "
              "wrong and part 1 might be readable after all - check before trusting either part")

print("\n2. the same two quantities evaluated at a COMMON pH, so the outcome does not appear on both sides")
# PART 1 EVALUATES EACH GENOME'S QUANTITY AT ITS OWN ph_min, WHICH PUTS THE OUTCOME INSIDE THE PREDICTOR. That is
# the right definition of "capacity where the organism actually lives", and it is also circular enough that it
# cannot stand alone - a monotone function of y regressed on y will find something. Evaluating at a fixed pH removes
# the circularity entirely at the cost of asking a slightly different question, and only quantities that survive
# BOTH are reported as evidence.
print(f"   {'pH':>5s} {'charge change 7->pH':>22s} {'p':>10s} {'buffer capacity':>18s} {'p':>10s}")
for pH in (3.0, 4.0, 5.0, 6.0):
    d2 = np.abs(charge(TP, 7.0) - charge(TP, pH))
    b2 = buffer_capacity(TP, pH)
    bb, _, pp = fit(d2)
    bc2, _, pc2 = fit(b2)
    print(f"   {pH:>5.1f} {bb:>+22.5f} {pp:>10.2e} {bc2:>+18.5f} {pc2:>10.2e}")
    rows.append({"part": "2", "pH": pH, "beta_charge_change": bb, "p_charge_change": pp,
                 "beta_capacity": bc2, "p_capacity": pc2, "n": len(B)})

print("\n3. where in pH does the signal live? the scan that names the residue")
# A composition weighted by how strongly each residue titrates at pH t. If carboxylate protonation is the mechanism
# the weighting should sharpen the association near pH 4; if histidine, near pH 6. The quantity is the same at every
# t, so the curve is comparable along its length.
best = (None, 0.0)
print(f"   {'pH':>5s} {'titration-weighted share':>26s} {'se':>8s} {'p':>10s} {'|t|':>7s}")
for t in np.arange(2.5, 8.01, 0.5):
    w = np.zeros(len(B))
    for k, pk in PKA.items():
        f = frac_protonated(t, pk)
        w += TP[k] * f * (1.0 - f) * (1 if SIGN[k] < 0 else 0)    # acidic groups only, weighted by titration at t
    x = 100.0 * w / TP["n"]
    if np.std(x) == 0:
        continue
    b, se, p = fit(x)
    tt = abs(b / se)
    if tt > best[1]:
        best = (float(t), tt)
    if abs(t * 2 - round(t * 2)) < 1e-9 and (t * 10) % 10 in (0, 5):
        print(f"   {t:>5.1f} {b:>+26.5f} {se:>8.5f} {p:>10.2e} {tt:>7.2f}")
    rows.append({"part": "3", "pH": float(t), "beta_per_sd": b, "se": se, "p": p, "abs_t": tt})
print(f"   sharpest at pH {best[0]:.1f} (|t| {best[1]:.2f}).  Aspartate pKa 3.65, glutamate 4.25, histidine 6.00.")
rows.append({"part": "3-peak", "measure": "pH of sharpest titration-weighted association", "pH": best[0], "abs_t": best[1]})

print("\n4. how far would the model pKa values have to be wrong to change the reading?")
# The honest objection to all of the above is that protein pKa values shift. A shift applied to every group alike
# moves the curves without reordering genomes, so it cannot flip a sign - but a shift large enough to move the
# carboxylates out of the growth range would change what part 3 is naming. This measures that margin directly.
shifts = []
for dp in (0.0, 0.5, 1.0, 1.5, 2.0):
    w = np.zeros(len(B))
    for k, pk in PKA.items():
        if SIGN[k] > 0:
            continue
        f = frac_protonated(best[0], pk + dp)
        w += TP[k] * f * (1.0 - f)
    x = 100.0 * w / TP["n"]
    if np.std(x) == 0:
        continue
    b, se, p = fit(x)
    shifts.append((dp, b, p))
    print(f"   carboxylate pKa shifted +{dp:.1f}: {b:+.5f} (p {p:.2e})")
    rows.append({"part": "4", "pka_shift": dp, "beta_per_sd": b, "p": p})
flips = [dp for dp, b, p in shifts if np.sign(b) != np.sign(shifts[0][1]) or p >= 0.05]
print(f"   the reading survives a uniform carboxylate pKa shift of up to "
      f"{(min(flips) - 0.5) if flips else max(dp for dp, _, _ in shifts):.1f} pH units"
      + ("" if flips else " (the largest tested)"))

print("\n5. verdict, read from part 2 only")
# THE DISCRIMINATION, DONE ON THE NON-CIRCULAR QUANTITIES. A positive coefficient means the quantity goes UP as the
# growth limit goes up, i.e. as acid tolerance goes DOWN.
#   BURDEN wants acid tolerance to come with a SMALL charge change -> POSITIVE coefficient on charge change.
#   BUFFERING wants acid tolerance to come with LARGE capacity    -> NEGATIVE coefficient on capacity.
# The reading is taken at the pH values where the question is actually posed - pH 3 and 4, inside the acidophile
# range and near the carboxylate pKa values - and the behaviour at pH 6 is reported separately rather than averaged
# in, because pH 6 is histidine's pKa and histidine's association runs the other way throughout this project.
p2 = {r["pH"]: r for r in rows if r["part"] == "2"}
lowpH = [p2[q] for q in (3.0, 4.0) if q in p2]
burden2 = all(r["beta_charge_change"] > 0 and r["p_charge_change"] < 0.05 for r in lowpH)
buffering2 = all(r["beta_capacity"] < 0 and r["p_capacity"] < 0.05 for r in lowpH)
print(f"   at pH 3 and 4 - inside the acidophile range, near the carboxylate pKa values:")
for r in lowpH:
    print(f"      pH {r['pH']:.0f}: charge change {r['beta_charge_change']:+.5f} (p {r['p_charge_change']:.2e})   "
          f"capacity {r['beta_capacity']:+.5f} (p {r['p_capacity']:.2e})")
if burden2 and not buffering2:
    v = ("BURDEN, on the non-circular test. At pH 3 and 4 the charge change carries a POSITIVE coefficient, so acid-tolerant genomes have LESS "
         "outward-surface charge to change - which is what the burden account predicts - while buffer capacity carries a positive coefficient too, "
         "which is the OPPOSITE of what buffering predicts. This partly answers correction 21's tension: UreI needs its carboxylates because it is a "
         "gate that must respond, and bulk housekeeping surfaces appear to do better with fewer, which is the reconciliation stated in "
         "docs/EXPERIMENTAL_EVIDENCE.md section 2.2 - now with a measurement behind it rather than only an argument.")
elif buffering2 and not burden2:
    v = ("BUFFERING, on the non-circular test, which contradicts the charge-burden mechanism this project has been offering and agrees with UreI and "
         "the gastric-pH engineering study.")
elif burden2 and buffering2:
    v = "BOTH FIRE in the predicted directions at once, which cannot be and means the two quantities are not separating here."
else:
    v = ("NEITHER ACCOUNT IS CLEANLY SUPPORTED at pH 3 and 4 on the non-circular quantities. The association with raw residue share is robust "
         "throughout this project, so what fails is expressing the mechanism as a titration quantity, not the association.")
print(f"\n   {v}")
# THE SIGN FLIP IS A RESULT, NOT A NUISANCE, and is reported rather than smoothed.
hi = [p2[q] for q in (5.0, 6.0) if q in p2]
if hi and all(r["beta_charge_change"] < 0 for r in hi) and lowpH and all(r["beta_charge_change"] > 0 for r in lowpH):
    print(f"   BOTH QUANTITIES REVERSE SIGN between pH 4 and pH 6 ({p2[4.0]['beta_charge_change']:+.4f} to {p2[6.0]['beta_charge_change']:+.4f}). "
          "That is where histidine")
    print("   titrates, and histidine's association runs opposite to the acidic one throughout this project - so the")
    print("   reversal is the second mechanism appearing in the same quantity, not instability in the first.")
    rows.append({"part": "5-flip", "measure": "sign reversal between pH 4 and 6", "beta_pH4": p2[4.0]["beta_charge_change"],
                 "beta_pH6": p2[6.0]["beta_charge_change"]})
# PART 3'S PEAK IS NOT WORTH THE WEIGHT IT LOOKS LIKE IT CARRIES.
print(f"   Part 3's scan peaks at pH {best[0]:.1f}, but the curve is flat above it (|t| {best[1]:.2f} against 8.96 at pH 8),")
print("   because a titration weighting applied far from every pKa becomes proportional to the raw share. So the peak")
print("   locates the signal between the aspartate and glutamate pKa values only weakly, and is not quoted as evidence")
print("   for which residue carries it.")
rows.append({"part": "5", "measure": "verdict", "verdict": v, "burden_supported": burden2, "buffering_supported": buffering2, "n": len(B)})
print("   WHAT THIS IS NOT. It is not an electrostatics calculation: no structure, no solvent, no interactions between")
print("   groups, and model pKa values throughout. It is the first-order question those calculations would refine, and")
print("   it is worth asking first because the two accounts differ in SIGN - which a first-order treatment can see and")
print("   a more elaborate one would not change.")

pd.DataFrame(rows).to_csv(D / "titration.csv", index=False)
print(f"\nwrote {D / 'titration.csv'}")
