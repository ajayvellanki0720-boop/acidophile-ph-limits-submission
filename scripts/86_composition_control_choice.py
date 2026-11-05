"""Does the headline survive a DIFFERENT composition control, or only the one it was built with? -> DATA_DIR/composition_control_choice.csv and, with `wide`, DATA_DIR/composition_control_choice_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
THE CLAIM THIS TESTS IS THE CENTRAL ONE. Everything here rests on the structural panel's exposed-periplasmic acidic share predicting the acid limit BEYOND what whole-proteome
composition predicts, and "whole-proteome composition" has been operationalised throughout as exactly one covariate: genome-wide acidic excess. That choice was never
compared against an alternative. The genome-feature table carries nine other charge measures, and scripts/85 found one of them - mean membrane-protein isoelectric point -
predicting the limit in its own right at p 3.6e-06, which is the specific reason to ask.
ADDING THEM IS THE WRONG TEST, AND SAYING WHY IS HALF THE POINT. Every one of these measures is a near-duplicate of the control already in the model: mean isoelectric point
correlates -0.953 with acidic excess, the share of proteins below pI 5 correlates +0.960, and even the most distant of the pI measures is at -0.848. Entering them beside it
would inflate both standard errors and invite reading a null as absence of effect, which is the error scripts/82 had to withdraw a docstring claim over. So the test is
SUBSTITUTION: swap the composition control for each alternative in turn and ask whether the panel coefficient survives. A claim that holds only under the one control it was
built with is a claim about that control.
AND ONE WELL-CONDITIONED WAY TO CONTROL THEM ALL AT ONCE. The first principal component of the standardised charge measures is a single column that carries what they share,
so it can be entered BESIDE acidic excess without the collinearity that defeats entering them individually. Its loadings and variance share are printed, because a composite
nobody can interpret is not a control anybody should trust.
THE CRITERION IS FIXED BEFORE THE NUMBERS, as in scripts/84: a sign change or a crossing of p=0.05 is a qualitative failure and would mean the claim depends on the choice of
control. A move within the coefficient's own standard error is quantitative, and the largest such move is the number to quote.
SCALE. Per percentage point of pooled exposed-periplasmic share, from ph.pgls directly and never ph.pgls_scan - correction 15. The other charged-residue shares of the same
class are held fixed in every fit, which is correction 14 and the error this project has now made three times in two days.
Usage: python3 scripts/86_composition_control_choice.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
CHARGE = ["acidic_excess", "mean_pI", "median_pI", "frac_pI_gt7", "frac_pI_lt5", "mean_pI_membrane",
          "mean_pI_soluble", "frac_acidic", "frac_basic"]
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
GF = ["ncbi_assembly"] + CHARGE + ["mean_gravy", "gc"]
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[GF], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]

tot = {k: np.zeros(len(B)) for k in ("n", "D", "E", "K", "R")}
for ref in MEM:
    if f"{ref}|peri_exp|n" not in B:
        continue
    tot["n"] += B[f"{ref}|peri_exp|n"].fillna(0).to_numpy(float)
    for k in ("D", "E", "K", "R"):
        tot[k] += B[f"{ref}|peri_exp|{k}"].fillna(0).to_numpy(float)
ok = tot["n"] >= MINPOS
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
y = g.ph_min.to_numpy(float)
SH = {k: 100 * tot[k][ok] / tot["n"][ok] for k in ("D", "E", "K", "R")}
TARGET = SH["D"] + SH["E"]                     # the headline measure: exposed-periplasmic acidic share
OTHERS = [SH["K"], SH["R"]]                    # correction 14: the other charged shares stay in
XS = np.column_stack([np.ones(len(g))] + [g[k] for k in SRCC])
print(f"sample: {int(ok.sum())} of {len(B)} bacteria clear {MINPOS} pooled exposed-periplasmic positions")

print("\n1. how near-duplicate these controls are, which is why the test is substitution and not addition")
print(f"   {'control':>22} {'r with acidic excess':>21} {'sd':>9}")
for c in CHARGE[1:]:
    r = float(np.corrcoef(g[c].to_numpy(float), g.acidic_excess.to_numpy(float))[0, 1])
    print(f"   {c:>22} {r:>+21.3f} {g[c].std():>9.4f}")
    rows.append({"part": "collinearity", "control": c, "r_with_acidic_excess": r, "sd": float(g[c].std())})
print("   Entering any of these beside acidic excess would inflate both standard errors; mean hydropathy, at")
print(f"   r {float(np.corrcoef(g.mean_gravy, g.acidic_excess)[0, 1]):+.3f}, is the only nearly independent one and is used that way in part 3.")


def fit(ctl_cols, label):
    """The headline measure with `ctl_cols` as the composition control(s)."""
    X = np.column_stack([XS] + list(ctl_cols) + OTHERS + [TARGET])
    assert np.linalg.matrix_rank(X) == X.shape[1], f"{label}: design matrix is rank-deficient"
    lam = ph.pgls(y, np.column_stack([XS] + list(ctl_cols)), Cm)["lambda"]
    r = ph.pgls(y, X, Cm, lam=lam)
    return float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1]), lam


base = fit([g.acidic_excess.to_numpy(float) * 100], "baseline")
# TIED TO scripts/84 RATHER THAN LEFT TO LOOK LIKE A THIRD NUMBER. This script holds the lysine and arginine
# shares fixed (correction 14), which scripts/84's composite rows do not, so the baseline here is that script's
# "D+E, holding K and R fixed" row and NOT its published +0.10681. The two differ by 1%, and the discrepancy is
# what led to finding that scripts/84's stated reason for leaving composites uncontrolled was wrong - and that
# holding histidine fixed as well costs the headline 27%, which is now reported in that script's part 9. An
# assertion keeps the two in step: three scripts quoting three figures for one measure is correction 12.
assert abs(base[0] - 0.10540) < 1e-5, (f"baseline {base[0]:+.5f} is not scripts/84's K-and-R-controlled D+E "
                                       "figure; one of the two scripts has changed specification")
print(f"\n2. SUBSTITUTION: one composition control at a time, the headline measure refitted under each")
print(f"   {'composition control':>22} {'beta/pt':>10} {'se':>9} {'p':>10} {'move/se':>8} {'lambda':>7}")
print(f"   {'acidic excess (as built)':>22} {base[0]:>+10.5f} {base[1]:>9.5f} {base[2]:>10.2e} {'-':>8} {base[3]:>7.3f}")
rows.append({"part": "substitution", "control": "acidic excess (as built)", "beta_per_pt": base[0],
             "se": base[1], "p": base[2], "lambda": base[3], "move_in_se": 0.0})
flips, worst, worstlab = [], 0.0, None
for c in CHARGE[1:]:
    b, se, pv, lam = fit([g[c].to_numpy(float)], c)
    mv = abs(b - base[0]) / base[1]
    flip = ((b > 0) != (base[0] > 0)) or ((pv < 0.05) != (base[2] < 0.05))
    if flip:
        flips.append((c, b, pv))
    if mv > worst:
        worst, worstlab = mv, c
    print(f"   {c:>22} {b:>+10.5f} {se:>9.5f} {pv:>10.2e} {mv:>8.2f} {lam:>7.3f}"
          + ("   FLIP" if flip else ""))
    rows.append({"part": "substitution", "control": c, "beta_per_pt": b, "se": se, "p": pv, "lambda": lam,
                 "move_in_se": mv, "flip": flip})

print("\n3. all of them at once, as a principal component, which is the only well-conditioned way to do it")
Z = np.column_stack([(g[c].to_numpy(float) - g[c].mean()) / g[c].std() for c in CHARGE])
U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
share = float(S[0] ** 2 / (S ** 2).sum())
pc1 = U[:, 0] * S[0]
print(f"   PC1 of the {len(CHARGE)} standardised charge measures carries {share:.1%} of their variance; loadings:")
print("      " + ", ".join(f"{c} {v:+.2f}" for c, v in zip(CHARGE, Vt[0])))
rows.append({"part": "pc1", "variance_share": share,
             "loadings": ";".join(f"{c}={v:+.3f}" for c, v in zip(CHARGE, Vt[0]))})
for lab, cols in (("PC1 alone", [pc1]),
                  ("PC1 + acidic excess", [pc1, g.acidic_excess.to_numpy(float) * 100]),
                  ("acidic excess + hydropathy", [g.acidic_excess.to_numpy(float) * 100, g.mean_gravy.to_numpy(float)]),
                  ("acidic excess + hydropathy + GC", [g.acidic_excess.to_numpy(float) * 100,
                                                       g.mean_gravy.to_numpy(float), g.gc.to_numpy(float)])):
    b, se, pv, lam = fit(cols, lab)
    mv = abs(b - base[0]) / base[1]
    flip = ((b > 0) != (base[0] > 0)) or ((pv < 0.05) != (base[2] < 0.05))
    if flip:
        flips.append((lab, b, pv))
    if mv > worst:
        worst, worstlab = mv, lab
    print(f"   {lab:>30} {b:>+10.5f} (se {se:.5f}, p {pv:.2e}), {mv:.2f} se from the baseline"
          + ("   FLIP" if flip else ""))
    rows.append({"part": "multiple controls", "control": lab, "beta_per_pt": b, "se": se, "p": pv,
                 "lambda": lam, "move_in_se": mv, "flip": flip})

print("\n4. the verdict, by the criterion fixed in the preamble")
n = len([r for r in rows if r["part"] in ("substitution", "multiple controls")])
print(f"   {n} specifications of the composition control; {len(flips)} change the headline measure's sign or")
print(f"   significance. Largest move {worst:.2f} standard errors ({worstlab}).")
for lab, b, pv in flips:
    print(f"      {lab}: {b:+.5f} (p {pv:.2e})")
rows.append({"part": "verdict", "n_specifications": n, "n_flips": len(flips),
             "largest_move_in_se": worst, "largest_move_control": worstlab})
assert not flips, ("the headline measure's sign or significance depends on which composition control is used. "
                   "That is not a caveat - it would mean the central claim is a claim about the control.")
print("\n5. the one specification that does move it, and it is not a composition control at all")
# The largest move in part 2 is mean membrane-protein pI at 1.61 standard errors, which is still within the
# criterion. But scripts/84 part 9 found a specification OUTSIDE this script's question that costs the headline
# 27%: holding exposed-periplasmic histidine fixed. That is worth placing beside these results, because the point
# of this script is that the claim does not depend on which composition control is chosen - and it does not - while
# the thing it DOES depend on is a residue share on the same surface, which is a different kind of dependence and
# a sharper question.
print("   Across 13 composition controls the headline moves at most 1.61 standard errors and never changes sign")
print("   or significance, so the central claim is not an artefact of the one control it was built with.")
print("   The specification that does cost it something is not a composition control: holding exposed-periplasmic")
print("   HISTIDINE fixed takes it from +0.10681 to +0.07745, a 27% reduction, still at p 5.9e-09 (scripts/84")
print("   part 9). Almost all of that is histidine rather than the basic residues, which cost 1% between them.")
print("   So the open question this leaves is not about whole-proteome composition - which is now well covered -")
print("   but about whether the exposed acidic share is partly standing in for exposed histidine on the same")
print("   surface, or whether the two are separate mechanisms that co-vary. Neither script can tell them apart.")
rows.append({"part": "pointer", "note": "holding exposed-periplasmic histidine fixed costs the headline 27% "
                                        "(+0.10681 -> +0.07745, p 5.9e-09); see scripts/84 part 9",
             "largest_composition_move_in_se": worst, "largest_composition_control": worstlab})

pd.DataFrame(rows).to_csv(D / f"composition_control_choice{TAG}.csv", index=False)
print("\nwrote", D / f"composition_control_choice{TAG}.csv")
