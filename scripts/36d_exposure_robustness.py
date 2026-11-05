"""Does the exposed-periplasmic signal survive the tests that killed earlier leads? -> DATA_DIR/exposure_robustness.csv
Every test uses one composite predictor: the mean of the z-scored exposed-periplasmic acidic excess over the membrane references a genome can be scored for (YidC, SecY from both
reference species, ATP synthase c, NuoL, NuoN, CyoB), which is more robust than any single family. The matched control predictors are built the same way from the buried
periplasmic, exposed cytoplasmic and transmembrane classes and from the soluble-reference surfaces, so every claim of specificity is like-for-like.
Tests: (1) PGLS on the bacterial marker tree, with and without genome-wide acidic excess, GC and genome size as covariates; (2) WITHIN FAMILIES - family-mean-centred pH_min on
family-mean-centred predictor, significance by permuting the predictor inside each family (this is the test the histidine and Trk/Ktr leads failed); (3) by data source, so the
hand-curated literature rows can be seen on their own; (4) by bound type, keeping only minima recorded as measured.
Usage: python3 scripts/36d_exposure_robustness.py [n_permutations]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
SOL = ["EF-Tu", "GroEL", "RecA"]
reflen = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
B["gcz"] = (B.gc - B.gc.mean()) / B.gc.std()
B["sz"] = (B.genome_size_bp - B.genome_size_bp.mean()) / B.genome_size_bp.std()


def composite(refs, cls):
    """mean of the z-scored acidic excess of `cls` over the references that a genome can be scored for"""
    cols = []
    for ref in refs:
        for c in ([cls] if isinstance(cls, str) else cls):
            if f"{ref}|{c}|n" not in B or (ref, c) not in reflen:
                continue
            n = B[f"{ref}|{c}|n"]
            x = ((B[f"{ref}|{c}|D"] + B[f"{ref}|{c}|E"] - B[f"{ref}|{c}|K"] - B[f"{ref}|{c}|R"]) / n).where(n >= max(8, 0.5 * reflen[(ref, c)]))
            cols.append((x - x.mean()) / x.std())
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3), M.notna().sum(axis=1)


PRED = {"exposed periplasmic": composite(MEM, "peri_exp")[0], "buried periplasmic": composite(MEM, "peri_bur")[0],
        "exposed cytoplasmic": composite(MEM, "cyto_exp")[0], "transmembrane": composite(MEM, "TM")[0], "soluble surface": composite(SOL, "surface")[0]}
rows = []


def pgls(mask, x, cov):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    use = [c for c in cov if g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, np.asarray(x, float)[:, None])
    return b[0], se[0], p[0], len(g)


SRCC = [f"src_{o}" for o in SRC]
for name, pr in PRED.items():
    for lab, cov in (("plain", SRCC), ("+ genome-wide acidic excess", ["gw"] + SRCC), ("+ acidic excess, GC, genome size", ["gw", "gcz", "sz"] + SRCC), ("+ optimum", ["ph_opt"] + SRCC)):
        m = pr.notna().to_numpy() & (B.ph_opt.notna().to_numpy() if lab == "+ optimum" else True)
        if m.sum() < 200:
            continue
        b, se, p, n = pgls(m, pr[m].to_numpy(), cov)
        rows.append({"test": "PGLS, all bacteria", "predictor": name, "model": lab, "n": n, "beta_per_sd": b, "se": se, "p": p})

# ---- within families ----
rng = np.random.default_rng(20261025)
for name, pr in PRED.items():
    d = B.assign(x=pr)[["family", "ph_min", "x", "gw"]].dropna()
    d = d[d.family.isin(d.family.value_counts()[lambda s: s >= 4].index)]
    if len(d) < 200:
        continue
    for lab, resid in (("raw", False), ("genome-wide acidic excess removed", True)):
        dd = d.copy()
        if resid:  # take out the between-genome composition effect first, within families
            A = np.column_stack([np.ones(len(dd)), dd.gw])
            dd["x"] = dd.x - A @ np.linalg.lstsq(A, dd.x, rcond=None)[0]
            dd["ph_min"] = dd.ph_min - A @ np.linalg.lstsq(A, dd.ph_min, rcond=None)[0]
        c = dd.groupby("family")[["ph_min", "x"]].transform(lambda s: s - s.mean())
        obs = np.polyfit(c.x, c.ph_min, 1)[0]
        null = np.empty(NPERM)
        gidx = [np.flatnonzero((dd.family == f).to_numpy()) for f in dd.family.unique()]
        xv, yv = c.x.to_numpy(), c.ph_min.to_numpy()
        for i in range(NPERM):
            xp = xv.copy()
            for ix in gidx:
                xp[ix] = rng.permutation(xp[ix])
            null[i] = np.polyfit(xp, yv, 1)[0]
        p = (1 + np.sum(np.abs(null) >= abs(obs))) / (NPERM + 1)
        rows.append({"test": f"within families ({dd.family.nunique()} families)", "predictor": name, "model": lab, "n": len(dd), "beta_per_sd": obs, "se": null.std(), "p": p})

# ---- by data source and bound type ----
pr = PRED["exposed periplasmic"]
for lab, m in (("hand-curated literature only", B.origin == "literature"), ("parser-extracted descriptions only", B.origin == "literature_auto"),
               ("databases only (ThermoBase/AciDB/BacDive)", B.origin.isin(["thermobase", "acidb", "bacdive"])), ("minima recorded as measured", B.ph_min_bound == "measured")):
    mm = (pr.notna() & m).to_numpy()
    if mm.sum() < 60:
        continue
    for sub, cov in (("plain", SRCC), ("+ genome-wide acidic excess", ["gw"] + SRCC)):
        b, se, p, n = pgls(mm, pr[mm].to_numpy(), cov)
        rows.append({"test": "subset: " + lab, "predictor": "exposed periplasmic", "model": sub, "n": n, "beta_per_sd": b, "se": se, "p": p})

R = pd.DataFrame(rows); R.to_csv(D / "exposure_robustness.csv", index=False)
pd.set_option("display.width", 200)
for t, g in R.groupby("test", sort=False):
    print(f"\n=== {t}   (beta per SD of the composite predictor)")
    print(g[["predictor", "model", "n", "beta_per_sd", "se", "p"]].round(4).to_string(index=False))
