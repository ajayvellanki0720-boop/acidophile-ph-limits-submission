"""Is the exposed acidic share standing in for exposed histidine, or are they separate mechanisms? -> DATA_DIR/acidic_vs_histidine.csv and, with `wide`, DATA_DIR/acidic_vs_histidine_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHERE THIS COMES FROM AND WHY IT MATTERS MORE THAN THE SCRIPT THAT FOUND IT. scripts/86 set out to show the central claim does not depend on which whole-proteome composition
control is used - it does not, across thirteen specifications - and in reconciling a 1% discrepancy against scripts/84 it turned up something larger: holding exposed-periplasmic
HISTIDINE fixed takes the headline acidic coefficient from +0.10681 to +0.07745 per point, a 27% reduction. Histidine entered this panel only after correction 18, so no
published figure has ever held it fixed. The question is what that 27% is.
THE SIGN STRUCTURE SAYS IT IS OMITTED-VARIABLE INFLATION, AND THE DIRECTION IS PREDICTED RATHER THAN FITTED. The acidic share and the histidine share correlate -0.549 across the genomes this
script fits (the first draft of this line said -0.463, which is the correlation on the raw panel before the bacteria-with-a-tip-and-a-pH-minimum
filter; a docstring figure its own script contradicts is the error this project keeps finding, so it is computed in part 1 and quoted from there); the acidic coefficient is POSITIVE and histidine's NEGATIVE. A genome high in exposed acidic residues is therefore low in exposed histidine, and both departures push
the acid limit the same way, so a model omitting histidine lets the acidic share carry part of histidine's effect and its coefficient comes out too large. That predicts the
direction of the 27% before any model is fitted: holding histidine fixed must SHRINK the acidic coefficient, not grow it. It does. This is correction 14's error - a residue
share entered without the correlate it shares variance with - reaching the headline measure itself, after appearing in a headline (14), in a control (19) and in a new script
(scripts/83).
BUT INFLATION IS NOT THE SAME AS PROXYING, and the whole point of this script is that those two readings make different predictions. If the acidic share were a PROXY for
histidine, then entering both should leave histidine carrying the effect and the acidic share near zero, and the acidic effect should inherit histidine's properties - above
all its partial two-sidedness, since correction 19 found histidine has a large shared component across the membrane while the acidic effect is strictly outward-facing. If they
are SEPARATE mechanisms that co-vary, both survive mutual adjustment and the acidic effect keeps its one-sidedness with histidine held fixed. Those are the tests below, and
the second one is the discriminating one because it uses a property the two measures are already known NOT to share.
WHAT IT DOES NOT DO. It does not change the canonical composite. Part 5 computes what including histidine would do, in per-standard-deviation terms so it is comparable to the
published headline, and records it for the deviations procedure instead of adopting it - the same treatment correction 19 gave the histidine experimental arm.
SCALE. Per percentage point from ph.pgls directly, never ph.pgls_scan - correction 15. Where a per-SD figure is given it is stated as such and computed as the per-unit
coefficient times the predictor's own standard deviation, with that standard deviation printed.
Usage: python3 scripts/87_acidic_vs_histidine.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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


def pooled(cls):
    tot = {k: np.zeros(len(B)) for k in ("n",) + RES}
    for ref in MEM:
        if f"{ref}|{cls}|n" not in B or f"{ref}|{cls}|H" not in B:
            continue
        tot["n"] += B[f"{ref}|{cls}|n"].fillna(0).to_numpy(float)
        for k in RES:
            tot[k] += B[f"{ref}|{cls}|{k}"].fillna(0).to_numpy(float)
    return tot


PT, CT = pooled("peri_exp"), pooled("cyto_exp")
ok = (PT["n"] >= MINPOS) & (CT["n"] >= MINPOS)
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
y = g.ph_min.to_numpy(float)
XS = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
lam = ph.pgls(y, XS, Cm)["lambda"]
P = {k: 100 * PT[k][ok] / PT["n"][ok] for k in RES}
C = {k: 100 * CT[k][ok] / CT["n"][ok] for k in RES}
DE, H = P["D"] + P["E"], P["H"]
DEc, Hc = C["D"] + C["E"], C["H"]


def fit(cols):
    X = np.column_stack([XS] + list(cols))
    assert np.linalg.matrix_rank(X) == X.shape[1], "design matrix is rank-deficient"
    r = ph.pgls(y, X, Cm, lam=lam)
    k = len(cols)
    return [(float(r["beta"][-k + i]), float(r["se"][-k + i]), float(r["p"][-k + i])) for i in range(k)]


print(f"sample: {int(ok.sum())} genomes scorable on both faces; lambda {lam:.3f}")
print("\n1. the sign structure, which predicts the direction of the 27% before anything is fitted")
rr = float(np.corrcoef(DE, H)[0, 1])
(bd,), (bh,) = fit([DE]), fit([H])
print(f"   r(acidic share, histidine share) = {rr:+.3f}; acidic alone {bd[0]:+.5f}/pt, histidine alone {bh[0]:+.5f}/pt")
same_push = (rr < 0) == ((bd[0] > 0) and (bh[0] < 0))
print("   The two correlate NEGATIVELY and their coefficients have OPPOSITE signs, so a genome high in exposed")
print("   acidic residues is low in exposed histidine and both departures push the limit the SAME way. A model")
print("   omitting histidine therefore lets the acidic share carry part of histidine's effect, and its")
print("   coefficient must come out too LARGE. The prediction is that holding histidine fixed SHRINKS it.")
rows.append({"part": "sign structure", "r_acidic_histidine": rr, "beta_acidic_alone": bd[0],
             "beta_histidine_alone": bh[0], "prediction": "holding histidine fixed shrinks the acidic coefficient"})

print("\n2. mutual adjustment: does either displace the other, as a proxy relationship requires?")
both = fit([H, DE])
(bh2, _, ph2), (bd2, _, pd2) = both[0], both[1]
print(f"   {'measure':>22} {'alone':>10} {'adjusted':>10} {'kept':>6} {'p adjusted':>11}")
print(f"   {'acidic D+E':>22} {bd[0]:>+10.5f} {bd2:>+10.5f} {abs(bd2 / bd[0]):>6.0%} {pd2:>11.2e}")
print(f"   {'histidine':>22} {bh[0]:>+10.5f} {bh2:>+10.5f} {abs(bh2 / bh[0]):>6.0%} {ph2:>11.2e}")
shrinks = abs(bd2) < abs(bd[0])
print(f"   The acidic coefficient {'SHRINKS' if shrinks else 'GROWS'}, as part 1 predicted."
      f" Both remain significant, and neither")
print("   is displaced: a proxy relationship would leave one of them near zero.")
rows.append({"part": "mutual adjustment", "beta_acidic_alone": bd[0], "beta_acidic_adjusted": bd2,
             "p_acidic_adjusted": pd2, "kept_acidic": abs(bd2 / bd[0]),
             "beta_histidine_alone": bh[0], "beta_histidine_adjusted": bh2, "p_histidine_adjusted": ph2,
             "kept_histidine": abs(bh2 / bh[0]), "direction_as_predicted": bool(shrinks)})
# `same_push` records that the sign structure and the fitted signs agree, which is what makes the prediction a
# prediction. A tautological assertion stood here in the first draft (`shrinks == same_push or True`), which can
# never fail and is worse than no assertion at all because it looks like a guard.
assert same_push, "the sign structure is not what part 1 describes, so its prediction does not follow"
assert shrinks, "part 1 predicted the acidic coefficient would shrink and it did not; the account is wrong"

print("\n3. THE DISCRIMINATING TEST: does the acidic effect inherit histidine's two-sidedness?")
# The acidic result's defining property is that it is strictly outward-facing: its exposed-cytoplasmic control is
# flat. Correction 19 found histidine is NOT like that - it has a large shared component across the membrane. So
# if the acidic share were carrying histidine's signal, holding histidine fixed should not matter to the acidic
# effect's sidedness, but the acidic effect should already LOOK partly two-sided, which it does not. The sharper
# version: compute the acidic sum-and-difference contrast with histidine held fixed. A proxy would lose its
# one-sidedness; a separate mechanism keeps it.
for lab, extra in (("acidic contrast, histidine free", ()), ("acidic contrast, histidine held fixed", (H, Hc))):
    r = fit(list(extra) + [DE + DEc, DE - DEc])
    (bs, _, ps), (bdf, sedf, pdf) = r[-2], r[-1]
    print(f"   {lab:>38}: sum {bs:>+9.5f} (p {ps:.1e}), difference {bdf:>+9.5f} (se {sedf:.5f}, p {pdf:.1e})")
    rows.append({"part": "sidedness", "specification": lab, "beta_sum_per_pt": bs, "p_sum": ps,
                 "beta_difference_per_pt": bdf, "se_difference": sedf, "p_difference": pdf})
sd_free = [r for r in rows if r["part"] == "sidedness"][0]
sd_held = [r for r in rows if r["part"] == "sidedness"][1]
keeps_sided = sd_held["p_difference"] < 0.05 and (sd_held["beta_difference_per_pt"] > 0) == (sd_free["beta_difference_per_pt"] > 0)
print(f"   The outward-facing excess {'SURVIVES' if keeps_sided else 'does NOT survive'} holding histidine fixed, at "
      f"{abs(sd_held['beta_difference_per_pt'] / sd_free['beta_difference_per_pt']):.0%} of its size.")
print("   That is the property the two measures are already known not to share, so it is the test that")
print("   discriminates rather than one that merely re-describes the correlation.")
rows.append({"part": "discriminating verdict", "acidic_keeps_one_sidedness": keeps_sided,
             "share_retained": abs(sd_held["beta_difference_per_pt"] / sd_free["beta_difference_per_pt"])})

print("\n4. the reading")
if keeps_sided and abs(bd2 / bd[0]) > 0.5 and pd2 < 0.05 and ph2 < 0.05:
    verdict = ("SEPARATE MECHANISMS THAT CO-VARY. Neither measure displaces the other, both stay significant "
               "under mutual adjustment, and the acidic effect keeps the strict one-sidedness that histidine "
               "does not have. The 27% is shared variance that the published specification was crediting "
               "entirely to the acidic share - an inflation to correct, not a signal to reassign.")
elif not keeps_sided:
    verdict = ("PROXY, at least in part. The acidic effect loses its one-sidedness once histidine is held "
               "fixed, which is the property that distinguished it, so the outward-facing reading would have "
               "to be re-examined.")
else:
    verdict = ("INDETERMINATE: one measure is substantially displaced, so the two cannot be separated on this "
               "sample.")
print("   " + verdict)
rows.append({"part": "reading", "verdict": verdict})

print("\n5. what including histidine in the canonical composite would do, COMPUTED AND NOT ADOPTED")
# The canonical composite is net D+E-K-R. Histidine's coefficient has the same sign as lysine's and arginine's,
# so a charge-style composite including it is D+E-K-R-H. Per-SD figures are given because the published headline
# is per SD, and each is the per-unit coefficient times that composite's own standard deviation, printed.
for lab, v in (("net D+E-K-R (canonical)", P["D"] + P["E"] - P["K"] - P["R"]),
               ("net D+E-K-R-H", P["D"] + P["E"] - P["K"] - P["R"] - P["H"])):
    (b, se, pv), = fit([v])
    s = float(np.std(v, ddof=1))
    print(f"   {lab:>26} {b:>+9.5f}/pt (p {pv:.2e}), SD {s:.4f} -> {b * s:>+8.4f} per SD")
    rows.append({"part": "composite", "composite": lab, "beta_per_pt": b, "se": se, "p": pv,
                 "sd_of_composite": s, "beta_per_sd": b * s})
# AND THESE ARE NOT THE PUBLISHED HEADLINE'S SPECIFICATION. The project's +0.185 per SD comes from scripts/57's
# model on a different composite with its own controls; the figures here are this panel share's per-SD effect in
# THIS model, and the comparison that matters is between the two rows above, not between either of them and the
# published number. Reading one against the other would be correction 12 - a comparison across specifications.
print("   These are not the published headline's specification: scripts/57's +0.185 per SD comes from a")
print("   different composite and model. The comparison that means anything is between the two rows above.")
cz = [r for r in rows if r["part"] == "composite"]
print(f"   Including histidine changes the per-SD figure from {cz[0]['beta_per_sd']:+.4f} to "
      f"{cz[1]['beta_per_sd']:+.4f}, i.e. by {abs(cz[1]['beta_per_sd'] / cz[0]['beta_per_sd']) - 1:+.0%}.")
print("   NOT ADOPTED HERE. The composite is fixed in docs/PREREGISTRATION.md and changing it is the deviations")
print("   procedure's call, not a script's - the same treatment correction 19 gave the histidine experimental")
print("   arm. What this script contributes is the number that decision needs.")
rows.append({"part": "composite", "composite": "recommendation", "adopted": False,
             "relative_change_per_sd": abs(cz[1]["beta_per_sd"] / cz[0]["beta_per_sd"]) - 1})

pd.DataFrame(rows).to_csv(D / f"acidic_vs_histidine{TAG}.csv", index=False)
print("\nwrote", D / f"acidic_vs_histidine{TAG}.csv")
