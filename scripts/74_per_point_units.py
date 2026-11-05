"""Every figure this project reports "per percentage point" is per STANDARD DEVIATION of that share -> DATA_DIR/per_point_units.csv, DATA_DIR/per_point_corrected.csv
CORRECTION 15, and the seventh instance of this project's one recurring error. `ph.pgls_scan` standardises the predictor it tests, so its coefficient is already per
standard deviation of the column as passed in - its own docstring says so in capitals, because four scripts had multiplied by the SD and converted it back to per unit.
What nobody checked is the opposite direction. SIX scripts pass a predictor measured in PERCENTAGE POINTS and store the returned coefficient in a column named
`beta_per_pt` or `pooled_beta_per_pt`: 34b, 34d, 34e, 36b, 36c and 40c. None divides by the standard deviation afterwards, so every one of those columns, and every
figure quoted from them, is per SD of that class's percentage share and not per percentage point. The factor is the share's own SD, which runs from about 2.5 to 6.5
points depending on class and panel, so the per-point effect is roughly three to six times SMALLER than the documents have been saying.
WHAT IS AND IS NOT WRONG, because the distinction decides how much has to be redone.
  - THE VALUES ARE CORRECT and the LABEL IS WRONG. Nothing was computed incorrectly; a per-SD coefficient was called per-point.
  - COMPARISONS BETWEEN CLASSES SURVIVE. The specificity table sets exposed periplasmic against cytoplasmic, buried, transmembrane and soluble surfaces, and per SD is a
    defensible common scale for that - it is the scale the headline composite uses. Those conclusions stand and only their units change.
  - ANY ARGUMENT ABOUT A PHYSICAL CHANGE AT A POSITION DOES NOT SURVIVE, and that is the damage. Flipping one of N positions changes a share by 1/N of 100 points, so
    reasoning from a position count to a predicted shift needs per-UNIT coefficients. scripts/73 did exactly that, differencing the D+E and lysine coefficients to
    compare a charge reversal against a neutralisation, and those two have different SDs, so it differenced two numbers on two scales. Its ratios are recomputed here.
THE CORRECTED COMPARISON STRENGTHENS CORRECTION 14 RATHER THAN WEAKENING IT. On the per-unit scale the margin reversal buys over neutralisation is 2.62x on the
hand-picked panel against 1.27x on the rule-selected one, a wider gap than the mislabelled figures showed, and the share of the reversal's shift that neutralisation
alone should recover is 38% against 79% - a 40-point gap in the charge-neutral control's expected position, not the 25 points first reported. Pooled per unit, lysine on
the rule-selected panel is -0.0024 at p 0.10, so the basic term is not merely weaker there but not significant at all.
Usage: python3 scripts/74_per_point_units.py   (env DATA_DIR=results/broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
rows = []
# Verified one by one rather than taken from the grep that first produced this list: the grep also returned
# scripts/29, which DOES divide (beta / sd_h) and is a second genuine exception alongside scripts/32b.
MISLABELLED = ["34b", "34d", "34e", "36b", "36c", "40c"]
EXCEPTIONS = {"32b": "uses ph.pgls directly, so its per-point column really is per point",
              "29": "divides the scan coefficient by the share SD before storing it"}


def _check_scan_is_per_sd(n=600, seed=3):
    """The mechanism, asserted on synthetic data so it fails here rather than in a document. A predictor measured
    on an arbitrary scale gets the SAME standardised coefficient from pgls_scan however it is rescaled, which is
    exactly what makes the stored column per SD rather than per unit."""
    rng = np.random.default_rng(seed)
    C = np.eye(n) + 0.25
    x = rng.normal(5.0, 4.0, n)
    y = 0.01 * x + rng.normal(0, 1.0, n)
    X = np.ones((n, 1))
    b1, *_ = ph.pgls_scan(C, y, X, x[:, None])
    b2, *_ = ph.pgls_scan(C, y, X, (x * 37.0)[:, None])      # rescaling the predictor must not move a per-SD beta
    assert abs(b1[0] - b2[0]) < 1e-9, f"pgls_scan is not scale-free: {b1[0]} vs {b2[0]}"
    direct = ph.pgls(y, np.column_stack([X, x]), C)["beta"][-1]
    assert abs(b1[0] / direct - x.std()) < 1e-6, "scan/direct should equal the predictor's SD exactly"
    return float(b1[0]), float(direct), float(x.std())


bs, bu, sd = _check_scan_is_per_sd()
print("0. the mechanism, checked on synthetic data before any real figure is touched")
print(f"   pgls_scan beta {bs:+.5f} | direct per-unit beta {bu:+.5f} | their ratio {bs / bu:.4f} = the predictor's SD {sd:.4f}")
print("   and multiplying the predictor by 37 leaves the scan's beta unchanged, which is what per-SD means.")
rows.append({"part": "mechanism", "scan_beta": bs, "direct_per_unit_beta": bu, "ratio": bs / bu, "predictor_sd": sd})


def pool(b, s):
    b, s = np.asarray(b, float), np.asarray(s, float)
    w = 1 / s ** 2
    fe = (w * b).sum() / w.sum()
    Q = (w * (b - fe) ** 2).sum()
    k = len(b) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (s ** 2 + tau2)
    mu = (wr * b).sum() / wr.sum()
    sem = float(np.sqrt(1 / wr.sum()))
    return float(mu), sem, float(2 * stats.norm.sf(abs(mu / sem)))


print(f"\n1. which stored columns carry the wrong label ({len(MISLABELLED)} scripts)")
print("   Every one computes with pgls_scan on a predictor in percentage points and stores the result unchanged.")
for s in MISLABELLED:
    rows.append({"part": "mislabelled column", "script": s,
                 "note": "stores a per-SD-of-share coefficient in a column named per_pt"})
print("   " + ", ".join(f"scripts/{s}" for s in MISLABELLED))
for s, why in EXCEPTIONS.items():
    print(f"   scripts/{s} is NOT affected: {why}.")
    rows.append({"part": "exception", "script": s, "note": why})
print("   The list was first recorded as seven and included scripts/29, which the grep matched because it uses")
print("   pgls_scan and stores a per-point column. Checking each by hand removed it: it divides. An unverified count")
print("   inside a correction is the kind of thing this project finds three days later.")

print("\n2. the residue decomposition on both scales, pooled across references by random effects")
PANELS = {"7-reference hand-picked": "", "32-reference rule-selected": "_wide"}
COMPS = ["D", "E", "D+E", "K", "R", "K+R", "net D+E-K-R"]
per_unit = {}
for pname, tag in PANELS.items():
    f = D / f"exposure_decomposition{tag}.csv"
    if not f.exists():
        print(f"   {pname}: {f.name} missing, skipped")
        continue
    d = pd.read_csv(f)
    # per unit = per SD divided by the share's own SD, which scripts/36c already stores per reference and component
    d["b_unit"] = d.beta_per_pt / d.sd_pts
    d["se_unit"] = d.se / d.sd_pts
    print(f"   {pname}")
    print(f"      {'component':14s} {'as published (per SD)':>22s} {'per POINT (correct)':>22s} {'mean share SD':>14s}")
    pu = {}
    for comp in COMPS:
        y = d[d.component == comp]
        if len(y) < 2:
            continue
        m_sd, s_sd, p_sd = pool(y.beta_per_pt, y.se)
        m_u, s_u, p_u = pool(y.b_unit, y.se_unit)
        pu[comp] = m_u
        print(f"      {comp:14s} {f'{m_sd:+.4f}':>22s} {f'{m_u:+.5f} (p {p_u:.0e})':>22s} {y.sd_pts.mean():>14.2f}")
        rows.append({"part": "both scales", "panel": pname, "component": comp, "n_refs": int(len(y)),
                     "pooled_beta_per_sd_of_share": m_sd, "se_per_sd": s_sd, "p_per_sd": p_sd,
                     "pooled_beta_per_point": m_u, "se_per_point": s_u, "p_per_point": p_u,
                     "mean_share_sd_points": float(y.sd_pts.mean())})
    per_unit[pname] = pu
print("   The two columns differ by the share's SD, so no conclusion about which component is larger changes within")
print("   a panel. What changes is every statement of the form 'X per percentage point', and anything derived from")
print("   one, because the per-point values are three to six times smaller.")

print("\n3. scripts/73's ratios, recomputed on the per-unit scale they always needed")
print("   Flipping one of N positions moves the D+E share by -100/N points and the K share by +100/N, so the")
print("   predicted shift is (100/N)(b_K - b_DE) for a reversal and (100/N)(-b_DE) for a neutralisation, with b in")
print("   PER-POINT units. N cancels; the SDs do not, which is why the mislabelled version was wrong.")
print(f"   {'panel':30s} {'ratio as published':>20s} {'ratio corrected':>18s} {'neutral share':>15s}")
for pname, pu in per_unit.items():
    if not {"D+E", "K"} <= set(pu):
        continue
    b_de, b_k = pu["D+E"], pu["K"]
    rev, neu = b_k - b_de, -b_de
    ratio, share = rev / neu, neu / rev
    old = pd.read_csv(D / "reversal_basis.csv") if (D / "reversal_basis.csv").exists() else None
    prev = np.nan
    if old is not None:
        q = old[(old.part == "reversal vs neutralisation") & (old.panel == pname)]
        if len(q):
            prev = float(q.iloc[0].ratio_reversal_over_neutralise)
    print(f"   {pname:30s} {f'{prev:.2f}x' if prev == prev else 'n/a':>20s} {f'{ratio:.2f}x':>18s} {share:>15.0%}")
    rows.append({"part": "corrected ratio", "panel": pname, "beta_DE_per_point": b_de, "beta_K_per_point": b_k,
                 "ratio_reversal_over_neutralise": ratio, "share_neutral_control_of_reversal": share,
                 "ratio_as_published_on_wrong_scale": prev})
if len(per_unit) == 2:
    a, b = (per_unit[k] for k in PANELS)
    s7 = -a["D+E"] / (a["K"] - a["D+E"])
    s32 = -b["D+E"] / (b["K"] - b["D+E"])
    print(f"   The charge-neutral control's expected position is {s7:.0%} against {s32:.0%}, a gap of "
          f"{abs(s32 - s7) * 100:.0f} percentage points rather than the 25 first reported.")
    print("   So correction 14 is strengthened by fixing the units, not weakened: the panels disagree more sharply")
    print("   about the control than the mislabelled figures suggested, and the hand-picked panel's rationale for")
    print("   reversing was correspondingly stronger than what replaced it.")
    rows.append({"part": "corrected neutral control", "share_7ref": s7, "share_32ref": s32,
                 "gap_percentage_points": abs(s32 - s7) * 100})

print("\n4. does the experiment's committed effect size depend on the mislabelled figures? No - but one claim does")
print("   docs/CHARGE_REVERSAL_DESIGN.md section 2 derives the predicted effect through the COMPOSITE, in SD units:")
print("   one reversal moves a protein's share by 2/n_positions, that is divided by the share's own SD to give a")
print("   z-shift, divided again by the number of references averaged, and multiplied by the per-SD composite slope.")
print("   That chain converts points to SDs explicitly and never multiplies by a per-point coefficient, so the")
print("   committed -0.05, the -0.16 optimistic end, and the n>=64 derived from them are UNAFFECTED by correction 15.")
# The exception is the row that row exists to reject: the naive per-reference transfer, which really does multiply a
# share change in points by the per-class pooled coefficient - mislabelled, so inflated by that share's SD.
NAIVE = {"NuoN": (4.65, 5), "NuoL": (2.35, 9)}      # points moved per reversal, and reversals in the primary arm
b_sd = next((r["pooled_beta_per_sd_of_share"] for r in rows
             if r.get("part") == "both scales" and r.get("panel", "").startswith("32") and r.get("component") == "net D+E-K-R"), np.nan)
b_pt = next((r["pooled_beta_per_point"] for r in rows
             if r.get("part") == "both scales" and r.get("panel", "").startswith("32") and r.get("component") == "net D+E-K-R"), np.nan)
print(f"\n   The naive row used the per-class pooled coefficient {b_sd:+.4f} as though it were per point. Per point it is")
print(f"   {b_pt:+.5f}, so that row is inflated by a factor of {b_sd / b_pt:.1f}:")
print(f"   {'protein':8s} {'pts/reversal':>13s} {'reversals':>10s} {'as published':>13s} {'corrected':>11s}")
tot_old = tot_new = 0.0
for prot, (pts, k) in NAIVE.items():
    old, new = pts * b_sd * k, pts * b_pt * k
    tot_old += old
    tot_new += new
    print(f"   {prot:8s} {pts:>13.2f} {k:>10d} {-old:>13.2f} {-new:>11.3f}")
    rows.append({"part": "naive transfer", "protein": prot, "points_per_reversal": pts, "n_reversals": k,
                 "naive_as_published": -old, "naive_corrected": -new})
print(f"   {'TOTAL':8s} {'':>13s} {sum(k for _, k in NAIVE.values()):>10d} {-tot_old:>13.2f} {-tot_new:>11.3f}")
print(f"   So the naive per-protein transfer predicts about {-tot_new:.2f} pH units for the primary arm's 14 reversals")
print(f"   against the composite basis's -0.05, so the gap that justifies using a composite is roughly")
print(f"   {tot_new / 0.05:.0f}-fold rather than the {tot_old / 0.05:.0f}-fold this reconstruction gives on the mislabelled scale.")
print(f"   This is a RECONSTRUCTION and not the design's own arithmetic: it reaches {-tot_old:.2f} where the document states")
print("   -1.9, because the document used each protein's own coefficient (NuoN +0.088) and this uses the pooled one.")
print("   The correction factor is the same either way - the share's SD - so the document's -1.9 becomes roughly -0.3.")
print("   The design's ARGUMENT survives - one protein of thousands should not deliver a whole coordinated shift, and")
print("   the composite basis is still the defensible one - but the word 'absurd' was carrying a mislabelled factor.")
rows.append({"part": "naive transfer total", "n_reversals": int(sum(k for _, k in NAIVE.values())),
             "naive_as_published": -tot_old, "naive_corrected": -tot_new, "composite_basis": -0.05,
             "gap_fold_as_published": tot_old / 0.05, "gap_fold_corrected": tot_new / 0.05,
             "note": "committed effect size is derived through the composite in SD units and is unaffected"})

R = pd.DataFrame(rows)
R.to_csv(D / "per_point_units.csv", index=False)
R[R.part == "both scales"].to_csv(D / "per_point_corrected.csv", index=False)
print("\nwrote", D / "per_point_units.csv", "and", D / "per_point_corrected.csv")
print("\nReading. Seventh instance of one error, and the first in the opposite direction: the four earlier ones multiplied")
print("a per-SD coefficient by an SD and reported per unit, this one left a per-SD coefficient labelled per unit. Both come")
print("from the same place, a wrapper that standardises silently, and the remedy is the same as correction 9's - the figure")
print("has to be derived a second way. It was found only because a ratio built on those coefficients had to be reasoned")
print("about physically, which forced the question of what one point of a share actually means. Conclusions about which")
print("class or residue matters are unaffected, because per SD is a fair common scale; every quoted PER-POINT number is not.")
