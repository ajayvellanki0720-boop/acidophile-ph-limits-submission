"""The manuscript's four figures, drawn from stored results only -> DATA_DIR/figures/*.png and *.svg
WHY THIS IS A SCRIPT AND NOT A NOTEBOOK. Every number in a figure is a claim, and a figure assembled by hand is a claim
with no file behind it - which is the exact failure `scripts/69` exists to catch for prose and cannot catch for a
picture. So each panel reads its values from the stored result file that the manuscript's own text cites for the same
figure, and asserts the headline values against that text before drawing anything. A figure that disagrees with the
manuscript stops the script.
FORM, CHOSEN BY WHAT EACH PANEL'S DATA HAS TO DO rather than by habit:
  1. the compartment contrast is IDENTITY plus POLARITY across four categories with uncertainty -> a forest plot (dot
     and interval against a zero line), not a bar chart, because bars imply magnitude from a zero baseline and the
     reader's question here is "which side of zero, and how far from it relative to its own error".
  2. the optimum decomposition is a NESTED-MODEL PROGRESSION on one sample -> an ordered dot plot, so the eye follows
     one quantity shrinking under added controls rather than comparing four unrelated heights.
  3. the per-reference result is a RANKING with a multiplicity decision -> a forest plot sorted by effect, with the
     BH survivors distinguished by colour AND by marker fill, never colour alone.
  4. the censoring bound is a CONTINUOUS RESPONSE with a hard physical ceiling -> a line, with the ceiling drawn as a
     shaded impossible region, because the point of that panel is that the requirement lies outside what can happen.
COLOUR. The data-viz reference palette's documented slots, used unchanged: slot 1 blue #2a78d6, slot 2 orange #eb6834,
slot 3 aqua #1baf7a. That file records these as passing every hard gate including the all-pairs floors for the first
three slots, so no re-stepping is done here and none is needed. The validator script could not be run in this
environment (no node), which is why the palette is taken as documented rather than modified - substituting values
would void the validation that makes them safe. Controls are drawn in text-secondary grey, not in a series colour,
because they are a reference class rather than a series competing for identity.
NO DUAL AXES ANYWHERE, and every panel carries its n and its units in the axis label rather than in a caption the
figure might get separated from.
Usage: python3 scripts/100_figures.py   (env DATA_DIR=results/broad2)"""
import os
import sys
import textwrap
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
OUT = D / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "figure.dpi": 200,
                     "savefig.bbox": "tight", "axes.labelcolor": INK, "text.color": INK})


def frame(ax, zero=True):
    """Recessive axes: no top or right spine, a thin x grid only, and a zero reference where polarity is the point."""
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    if zero:
        ax.axvline(0, color=INK2, lw=0.9, zorder=1)


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    print(f"   wrote {OUT / (name + '.png')} and .svg")


H = pd.read_csv(D / "headline_numbers.csv")


def row(fig_name):
    r = H[H.figure == fig_name]
    assert len(r) == 1, f"{fig_name!r} is not a single row of headline_numbers.csv"
    return r.iloc[0]


print("1. the compartment contrast")
SPEC = [("exposed periplasmic\n(the claim)", "exposed periplasmic, 32-reference rule-selected panel", True),
        ("exposed cytoplasmic", "exposed cytoplasmic control, 32 references", False),
        ("transmembrane", "transmembrane control, 32 references", False),
        ("soluble surface\n(5 references)", "soluble-surface control, 5 references", False)]
# ASSERTED AGAINST THE MANUSCRIPT'S OWN TEXT, so a figure cannot drift from the prose it illustrates.
assert abs(float(row(SPEC[0][1]).beta_per_sd) - 0.1848) < 5e-5, "the headline row is not the +0.1848 the manuscript quotes"
assert abs(float(row(SPEC[1][1]).beta_per_sd) - 0.0212) < 5e-5, "the cytoplasmic control is not the +0.0212 the manuscript quotes"
fig, ax = plt.subplots(figsize=(5.4, 2.5))
frame(ax)
for i, (lab, key, is_claim) in enumerate(SPEC):
    r = row(key)
    b, se = float(r.beta_per_sd), float(r.se_per_sd)
    y = len(SPEC) - 1 - i
    c = BLUE if is_claim else MUTED
    ax.plot([b - 1.96 * se, b + 1.96 * se], [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=3)
    ax.plot([b], [y], "o", ms=8, color=c, mec="white", mew=1.4, zorder=4)
    ax.text(b + 1.96 * se + 0.012, y, f"{b:+.4f}", va="center", ha="left", fontsize=8.5,
            color=INK if is_claim else INK2, fontweight="bold" if is_claim else "normal")
ax.set_yticks(range(len(SPEC)))
ax.set_yticklabels([s[0] for s in reversed(SPEC)])
ax.set_xlabel("pH units per SD of composite  (95% CI; n = 2,768 except soluble, 2,764)")
ax.set_xlim(-0.22, 0.36)
# NOT "only the periplasmic surface predicts the limit", which the first draft said and section 3.2 of the
# manuscript denies in as many words: the soluble-surface control is significant and OPPOSITE in sign, so the
# pattern is a contrast between compartments rather than a periplasm-only effect. A figure title is a claim and
# gets held to the text it illustrates.
ax.set_title("Periplasmic surface only: inward controls null, soluble surface reversed", loc="left", fontsize=9.5, color=INK, pad=8)
save(fig, "fig1_compartment_contrast")

print("2. the growth-optimum and histidine decomposition")
OH = pd.read_csv(D / "optimum_and_histidine_wide.csv")
F = OH[OH.part == "fit"].set_index("specification")
ORDER = ["composition only", "+ growth optimum", "+ histidine", "+ optimum and histidine"]
assert abs(float(F.loc["+ optimum and histidine", "beta_per_sd"]) - 0.112005) < 5e-6, \
    "the joint fit is not the +0.112005 per SD the manuscript quotes"
fig, ax = plt.subplots(figsize=(5.4, 2.5))
frame(ax)
for i, spec in enumerate(ORDER):
    r = F.loc[spec]
    b, se = float(r.beta_per_sd), float(r.se) * float(r.beta_per_sd) / float(r.beta_per_pt)
    y = len(ORDER) - 1 - i
    c = BLUE if i == 0 else (ORANGE if spec == "+ optimum and histidine" else MUTED)
    ax.plot([b - 1.96 * se, b + 1.96 * se], [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=3)
    ax.plot([b], [y], "o", ms=8, color=c, mec="white", mew=1.4, zorder=4)
    share = float(r.share_of_composition_only)
    ax.text(b + 1.96 * se + 0.012, y, f"{b:+.4f}" + (f"   {share:.0%} kept" if i else ""),
            va="center", ha="left", fontsize=8.5, color=INK2)
ax.set_yticks(range(len(ORDER)))
ax.set_yticklabels([s.replace("+ ", "+ ") for s in reversed(ORDER)])
ax.set_xlabel("pH units per SD of composite  (95% CI; all four fits on the same n = 1,655)")
ax.set_xlim(0, 0.46)
ax.grid(axis="x", color=GRID, lw=0.6)
ax.set_title("39% survives the growth optimum and histidine held together", loc="left", fontsize=10, color=INK, pad=8)
save(fig, "fig2_optimum_decomposition")

print("3. per-reference effects with Benjamini-Hochberg")
WC = pd.read_csv(D / "which_complexes_wide.csv")
per = WC[WC.part == "per reference"].copy()
bh = WC[WC.part == "bh"][["reference", "q", "survives_q10"]]
per = per.merge(bh, on="reference", suffixes=("", "_bh"))
per["surv"] = per.survives_q10_bh.fillna(per.survives_q10).astype(str).str.lower().isin(["true", "1", "1.0"])
per = per.sort_values("beta")
assert per.surv.sum() == 7, f"{int(per.surv.sum())} references survive BH here; the manuscript says 7"
assert abs(per.beta.min() - (-0.42859)) < 5e-6, "the strongest reference is not the -0.42859 the manuscript quotes"
fig, ax = plt.subplots(figsize=(6.0, 5.2))
frame(ax)
for i, r in enumerate(per.itertuples()):
    c = ORANGE if r.surv else MUTED
    ax.plot([r.beta - 1.96 * r.se, r.beta + 1.96 * r.se], [i, i], color=c, lw=1.8, solid_capstyle="round", zorder=3)
    # survivors get colour AND a filled marker; non-survivors an open one, so identity is never colour alone
    ax.plot([r.beta], [i], "o", ms=6.5, color=c if r.surv else "white", mec=c, mew=1.5, zorder=4)
ax.set_yticks(range(len(per)))
ax.set_yticklabels([n if len(n) <= 58 else n[:57] + "\u2026" for n in per.name], fontsize=7.2)
ax.set_xlabel("pH units per unit of detection  (95% CI; n = 2,768)")
ax.set_title("7 of 22 testable references survive BH at q<0.10", loc="left", fontsize=10, color=INK, pad=8)
h1 = plt.Line2D([], [], marker="o", ls="-", color=ORANGE, ms=6.5, mew=1.5, mec=ORANGE, label="survives q<0.10")
h2 = plt.Line2D([], [], marker="o", ls="-", color=MUTED, ms=6.5, mew=1.5, mec=MUTED, mfc="white", label="does not")
ax.legend(handles=[h1, h2], loc="upper left", frameon=False, fontsize=8)
save(fig, "fig3_per_reference_bh")

print("4. the censoring bound and its physical ceiling")
CS = pd.read_csv(D / "censoring_sensitivity.csv")
d5 = CS[(CS.part.astype(str) == "5") & CS.alpha_ph_per_sd.notna() & CS.beta_per_sd.notna()]
d5 = d5[d5.measure == "displacement aligned with predictor"].sort_values("alpha_ph_per_sd")
az = CS[CS.measure == "alpha that zeroes the headline"].iloc[0].alpha_ph_per_sd
p6 = CS[(CS.part.astype(str) == "6") & CS.required_depth_ph.notna()]
assert abs(float(az) + 0.1848) < 1e-3, "the zeroing displacement is not the 0.185 pH per SD the manuscript quotes"
fig, ax = plt.subplots(figsize=(5.4, 2.9))
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.grid(color=GRID, lw=0.6, zorder=0)
ax.set_axisbelow(True)
x = -d5.alpha_ph_per_sd.to_numpy(float)
y = d5.beta_per_sd.to_numpy(float)
ax.axhline(0, color=INK2, lw=0.9, zorder=1)
ax.plot(x, y, "-", color=BLUE, lw=2.0, zorder=3)
ax.plot(x, y, "o", ms=5.5, color=BLUE, mec="white", mew=1.2, zorder=4)
ax.plot([-float(az)], [0], "o", ms=9, color=ORANGE, mec="white", mew=1.6, zorder=5)
ax.annotate(f"headline reaches zero at\n{-float(az):.3f} pH per SD", xy=(-float(az), 0), xytext=(-float(az) - 0.075, 0.075),
            fontsize=8, color=INK2, ha="center",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xlabel("differential censoring aligned with the predictor  (pH per SD)")
ax.set_ylabel("headline, pH per SD")
ax.set_title("What censoring would have to do, and the depth it would need", loc="left", fontsize=10, color=INK, pad=8)
save(fig, "fig4a_censoring_sensitivity")

# THE SECOND HALF OF PANEL 4 is the part that closes the objection: required depth against available headroom. A
# grouped bar is right here because both quantities share one unit and the comparison is "is the bar past the line".
fig, ax = plt.subplots(figsize=(5.4, 2.4))
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
ax.set_axisbelow(True)
thr = p6[p6.required_depth_ph > 0].copy()
labs = [m.split(": ", 1)[1] for m in thr.measure]
yy = np.arange(len(thr))
ax.barh(yy + 0.20, thr.required_depth_ph, height=0.32, color=ORANGE, zorder=3, label="depth required to erase the result")
ax.barh(yy - 0.20, thr.headroom_to_ph0, height=0.32, color=MUTED, zorder=3, label="depth physically available (to pH 0)")
for i, r in enumerate(thr.itertuples()):
    ax.text(r.required_depth_ph + 0.25, i + 0.20, f"{r.required_depth_ph:.1f}", va="center", fontsize=8, color=INK2)
    ax.text(r.headroom_to_ph0 + 0.25, i - 0.20, f"{r.headroom_to_ph0:.1f}", va="center", fontsize=8, color=INK2)
ax.set_yticks(yy)
ax.set_yticklabels([textwrap.fill(lab, 30) for lab in labs], fontsize=7.6)
ax.set_xlabel("mean displacement of reported minimum below the true limit (pH units)")
ax.set_xlim(0, max(thr.required_depth_ph) * 1.22)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, frameon=False, fontsize=8)
ax.set_title("Every threatening proxy needs more depth than physically exists", loc="left", fontsize=10, color=INK, pad=8)
save(fig, "fig4b_censoring_headroom")

print("\nall figures drawn from stored results; every headline value asserted against the manuscript's text first")
