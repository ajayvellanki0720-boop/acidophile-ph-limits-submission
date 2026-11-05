"""Why do the within-family and between-species slopes disagree? -> DATA_DIR/slope_decomposition.csv and DATA_DIR/slope_decomposition_wide.csv
36d reported +0.418 pH per unit of the exposed-periplasmic composite within families (permutation test, polyfit on family-centred data) against +0.157 per SD from PGLS across
bacteria, and the sister-pair calibration in 38 favoured the smaller figure. The two numbers were never on the same scale: pgls_scan standardises the tested predictor
internally and so returns a slope per standard deviation, while the within-family figure came from polyfit and is a slope per UNIT of the composite, whose standard deviation is
not 1 (it is a mean of correlated z-scores) and which shrinks further once family means are removed. This script puts every estimator on both scales on one subset of genomes
and separates that scaling question from the genuine statistical ones: whether the within-family relationship differs from the between-family one (a Simpson-type pattern),
whether it survives phylogenetic correction inside families (family fixed effects in PGLS, the apples-to-apples comparison), and whether it rests on a few large families
(jackknife over families).
Usage: python3 scripts/39a_slope_decomposition.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
# panel selectable as in scripts/39b. The family-fixed-effects comparison here had only ever been run on the
# hand-picked panel, while the headline rests on the rule-selected one - the same gap scripts/39b had.
WHICH = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ("orig", "wide") else "orig"
if WHICH == "wide":
    _refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
    _refs = _refs[_refs.accepted & _refs.membrane]
    MEM = [str(r.label) for r in _refs.itertuples()]
    _classes, _panel = ROOT / "results/structure_panel_wide_classes.csv", D / "structure_panel_wide.csv"
else:
    MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
    _classes, _panel = ROOT / "results/structure_panel_classes.csv", D / "structure_panel.csv"
reflen = pd.read_csv(_classes).groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(_panel).drop_duplicates("ncbi_assembly"), on="ncbi_assembly"))
print(f"panel '{WHICH}': {len(MEM)} membrane references")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
cols = []
for r in MEM:
    if f"{r}|peri_exp|n" not in B:
        continue
    n = B[f"{r}|peri_exp|n"]
    x = ((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, "peri_exp")]))
    cols.append((x - x.mean()) / x.std())
M = pd.concat(cols, axis=1)
B["x"] = M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
B["gw"] = B.acidic_excess * 100
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
# one subset for every estimator: families with >= 4 species and both variables present
d = B[B.x.notna() & B.ph_min.notna()].copy()
d = d[d.family.isin(d.family.value_counts()[lambda s: s >= 4].index)].reset_index(drop=True)
keep = d.ncbi_assembly.tolist()
Call = ph.vcv_fast(tree, keep)
SRCC = [f"src_{o}" for o in SRC if d[f"src_{o}"].nunique() > 1]
print(f"{len(d)} genomes in {d.family.nunique()} families (>= 4 species each)\n")
sd_raw = d.x.std()
cx = d.x - d.groupby("family").x.transform("mean")
cy = d.ph_min - d.groupby("family").ph_min.transform("mean")
sd_cen = cx.std()
print("SCALE: the composite is a mean of correlated z-scores, so its SD is not 1")
print(f"  SD of composite across these genomes      {sd_raw:.3f}")
print(f"  SD of composite after removing family means {sd_cen:.3f}  ({sd_cen / sd_raw:.0%} of it)")
print(f"  share of variance that is BETWEEN families: composite {1 - (sd_cen / sd_raw) ** 2:.0%}, pH_min {1 - (cy.std() / d.ph_min.std()) ** 2:.0%}\n")
rows = []


def add(name, b_unit, se_unit, sdx, note=""):
    rows.append({"estimator": name, "n": len(d), "beta_per_unit": b_unit, "se_per_unit": se_unit, "sd_of_x_used": sdx,
                 "beta_per_sd_of_x_used": b_unit * sdx, "beta_per_sd_of_raw_composite": b_unit * sd_raw, "note": note})


def ols(y, x, extra=None):
    A = np.column_stack([np.ones(len(x)), np.asarray(x, float)] + ([np.asarray(e, float) for e in extra] if extra else []))
    bb, *_ = np.linalg.lstsq(A, np.asarray(y, float), rcond=None)
    resid = np.asarray(y, float) - A @ bb
    s2 = resid @ resid / (len(x) - A.shape[1])
    cov = s2 * np.linalg.pinv(A.T @ A)
    return bb[1], np.sqrt(cov[1, 1])


b, se = ols(d.ph_min, d.x); add("pooled OLS", b, se, sd_raw)
b, se = ols(d.ph_min, d.x, [d.gw]); add("pooled OLS + composition", b, se, sd_raw)
fm = d.groupby("family").agg(y=("ph_min", "mean"), x=("x", "mean"), g=("gw", "mean"))
b, se = ols(fm.y, fm.x); add("between families (family means)", b, se, fm.x.std(), f"{len(fm)} families")
b, se = ols(fm.y, fm.x, [fm.g]); add("between families + composition", b, se, fm.x.std(), f"{len(fm)} families")
b, se = ols(cy, cx); add("within families (OLS on centred data)", b, se, sd_cen, "the 36d estimator")
cg = d.gw - d.groupby("family").gw.transform("mean")
b, se = ols(cy, cx, [cg]); add("within families + composition", b, se, sd_cen, "the 36d conditioned estimator")
X0 = np.column_stack([np.ones(len(d))] + [d[c] for c in SRCC])
bb, ss, pp, lam = ph.pgls_scan(Call, d.ph_min.to_numpy(float), X0, d[["x"]].to_numpy(float))
add("PGLS (no family effects)", bb[0] / sd_raw, ss[0] / sd_raw, sd_raw, f"lambda {lam:.2f}")
Xg = np.column_stack([np.ones(len(d)), d.gw] + [d[c] for c in SRCC])
bb, ss, pp, lam = ph.pgls_scan(Call, d.ph_min.to_numpy(float), Xg, d[["x"]].to_numpy(float))
add("PGLS + composition", bb[0] / sd_raw, ss[0] / sd_raw, sd_raw, f"lambda {lam:.2f}")
dum = pd.get_dummies(d.family, drop_first=True).to_numpy(float)
Xf = np.column_stack([np.ones(len(d)), dum] + [d[c] for c in SRCC])
bb, ss, pp, lam = ph.pgls_scan(Call, d.ph_min.to_numpy(float), Xf, d[["x"]].to_numpy(float))
add("PGLS + family fixed effects", bb[0] / sd_raw, ss[0] / sd_raw, sd_cen, f"lambda {lam:.2f}; p {pp[0]:.1e}")
Xfg = np.column_stack([np.ones(len(d)), dum, d.gw] + [d[c] for c in SRCC])
bb, ss, pp, lam = ph.pgls_scan(Call, d.ph_min.to_numpy(float), Xfg, d[["x"]].to_numpy(float))
add("PGLS + family effects + composition", bb[0] / sd_raw, ss[0] / sd_raw, sd_cen, f"lambda {lam:.2f}; p {pp[0]:.1e}")
R = pd.DataFrame(rows)
R.insert(0, "panel", WHICH)
out = D / f"slope_decomposition{'_wide' if WHICH == 'wide' else ''}.csv"
R.to_csv(out, index=False)
pd.set_option("display.width", 230)
print("Every estimator on one subset. 'per unit' is the raw slope; the next two columns rescale it.")
print(R[["estimator", "beta_per_unit", "se_per_unit", "sd_of_x_used", "beta_per_sd_of_x_used", "beta_per_sd_of_raw_composite", "note"]].round(4).to_string(index=False))
# jackknife the within-family estimators over families
for nm, extra in (("within families (OLS)", None), ("within families + composition", [cg])):
    fams = d.family.unique(); est = []
    for f in fams:
        m = (d.family != f).to_numpy()
        yy = cy[m]; xx = cx[m]
        est.append(ols(yy, xx, [e[m] for e in extra] if extra else None)[0])
    est = np.array(est)
    k = np.argsort(np.abs(est - np.median(est)))[-3:]
    print(f"\njackknife over {len(fams)} families, {nm}: median {np.median(est):.3f}, range {est.min():.3f} to {est.max():.3f}; "
          f"most influential {', '.join(f'{fams[i]} ({est[i]:.3f})' for i in k)}")
