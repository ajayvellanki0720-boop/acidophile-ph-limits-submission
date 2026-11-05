"""Which references carry the signal once their own noise is removed, and is the experiment aimed at them? -> DATA_DIR/reference_ranking.csv
Script 53 found one dominant factor with loadings running from -0.09 to +0.43 against the +0.23 that equal contribution would give, which means the references are not
interchangeable. The charge-reversal design picked its target protein family before any of that was known, so the obvious question is whether it picked a good one.
Ranking references honestly needs one correction. A reference's own coefficient against the lower limit is attenuated by ITS OWN measurement noise, so a protein whose
alignment and structure happen to be cleaner will look more important than one that is equally relevant but noisier - which is exactly the confound that would make this
ranking useless for choosing an experimental target. Each reference's reliability is therefore estimated from its item-rest correlation: if every reference is the common
factor plus independent noise, corr(x_i, mean of the rest) = sqrt(rel_i * rel_rest), and rel_rest follows from the rest's mean inter-item correlation by Spearman-Brown. The
coefficient is then divided by sqrt(rel_i), which is the right power for a per-standard-deviation coefficient (see scripts/52 and correction 7 - this is the mistake that was
made once already and is simulated rather than asserted here too).
RESULT, recorded here because it determines how the output must be read: the disattenuation FAILS as a ranking device. Most references have very low individual reliability -
item-rest correlations put two thirds of them below 0.2 and several below 0.03 - so dividing by its square root divides by a number near zero and amplifies noise far more than
it corrects attenuation. A reference with a coefficient of -0.015 and a reliability of 0.002 comes out at -0.32, which is arithmetic, not evidence. So the disattenuated column
is reported but is only interpretable above a reliability floor, and the script applies one.
Two things do survive, and they are the point of the script:
  1. The RAW per-reference coefficient correlates about +0.8 with that reference's reliability. Any ranking of individual proteins by their own coefficient is therefore mostly
     a ranking of annotation and alignment quality, which is a warning against reading individual references' estimates as biology - including the ones this project has quoted.
  2. The family the charge-reversal design targets can be checked directly, and the answer is not the one the design assumed.
Usage: python3 scripts/54_reference_ranking.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.linalg import solve_triangular
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
NAME = {str(r.label): (str(r.name), str(r.pfam)) for r in refs.itertuples()}
DESIGN_PFAM = "PF00361"  # the NuoL/NuoM/NuoN family, which the charge-reversal design targets
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))

cols = {}
for r in MEM:
    c = f"{r}|peri_exp|n"
    if c not in B or (r, "peri_exp") not in reflen:
        continue
    n_ = B[c]
    x = ((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n_).where(n_ >= max(8, 0.5 * reflen[(r, "peri_exp")]))
    cols[r] = (x - x.mean()) / x.std()
Z = pd.DataFrame(cols)
names = list(Z.columns)
print(f"{len(names)} membrane references with a scorable exposed-periplasmic class\n")

# one factorisation per reference is unavoidable because each has its own scorable sample, but lambda is fitted once on the
# largest common subset and reused, since it describes the response's phylogenetic structure rather than the predictor's
base_keep = (Z.notna().sum(axis=1) >= 3) & B.ph_min.notna()
lam = ph.pgls(B.ph_min[base_keep].to_numpy(float),
              np.column_stack([np.ones(int(base_keep.sum()))] + [B.loc[base_keep, c].to_numpy(float) for c in ["gw"] + [f"src_{o}" for o in SRC]]),
              Vcv[np.ix_(base_keep.to_numpy(), base_keep.to_numpy())])["lambda"]
print(f"lambda fitted once on the {int(base_keep.sum())}-genome common subset: {lam:.3f}\n")


def fit_one(x):
    m = (x.notna() & B.ph_min.notna()).to_numpy()
    g = B[m]
    n = int(m.sum())
    C = Vcv[np.ix_(m, m)]
    Cn = C / np.mean(np.diag(C))
    L = np.linalg.cholesky(Cn * lam + np.diag(np.diag(Cn)) * (1 - lam) + 1e-10 * np.eye(n))
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X = np.column_stack([np.ones(n)] + [g[c].to_numpy(float) for c in use])
    Xs = solve_triangular(L, X, lower=True)
    Q, _ = np.linalg.qr(Xs)
    ys = solve_triangular(L, g.ph_min.to_numpy(float), lower=True)
    yt = ys - Q @ (Q.T @ ys)
    v = x[m].to_numpy(float)
    Fs = solve_triangular(L, (v - v.mean()) / v.std(), lower=True)
    Ft = Fs - Q @ (Q.T @ Fs)
    xx = float(Ft @ Ft)
    b = float(Ft @ yt) / xx
    s2 = (float(yt @ yt) - b * b * xx) / (n - Xs.shape[1] - 1)
    return b, np.sqrt(s2 / xx), n


load = pd.read_csv(D / "factor_structure.csv")
load = load[load.test == "loading"].set_index("ref").value.to_dict() if "ref" in load.columns else {}
rows = []
for r in names:
    b, se, n = fit_one(Z[r])
    rest = Z[[c for c in names if c != r]]
    pair = Z[r].corr(rest.mean(axis=1))
    tri = rest.corr().to_numpy()[np.triu_indices(len(names) - 1, 1)]
    rbar = float(np.nanmean(tri))
    k = len(names) - 1
    rel_rest = k * rbar / (1 + (k - 1) * rbar) if rbar > 0 else np.nan
    rel = min(1.0, pair ** 2 / rel_rest) if rel_rest and rel_rest > 0 and pair > 0 else np.nan
    dis = b / np.sqrt(rel) if rel and rel > 0 else np.nan
    nm, pf = NAME.get(r, ("?", "?"))
    rows.append({"ref": r, "name": nm, "pfam": pf, "n": n, "beta_per_sd": b, "se": se, "item_rest_r": float(pair),
                 "reliability": rel, "disattenuated": dis, "pc1_loading": load.get(r, np.nan)})

R = pd.DataFrame(rows).sort_values("disattenuated", ascending=False)
R.to_csv(D / "reference_ranking.csv", index=False)
print(f"reliability of the rest of the panel (Spearman-Brown on {len(names) - 1} references): {rel_rest:.3f}")
print(f"{'reference':9s} {'beta/SD':>8} {'rel':>6} {'disatt':>8} {'PC1':>7}  name")
for x in R.itertuples():
    star = "  <- the design's family" if x.pfam == DESIGN_PFAM else ""
    print(f"{x.ref:9s} {x.beta_per_sd:>+8.4f} {x.reliability if x.reliability == x.reliability else float('nan'):>6.3f} "
          f"{x.disattenuated if x.disattenuated == x.disattenuated else float('nan'):>+8.4f} "
          f"{x.pc1_loading if x.pc1_loading == x.pc1_loading else float('nan'):>+7.3f}  {x.name[:46]}{star}")

FLOOR = 0.20
ok = R[R.disattenuated.notna()]
rk = ok[ok.reliability >= FLOOR]
print(f"\n{len(ok)} references have an estimable reliability, but only {len(rk)} reach {FLOOR:.2f}.")
print(f"RAW coefficient against reliability: r = {ok[['beta_per_sd', 'reliability']].corr().iloc[0, 1]:+.3f}. A per-reference coefficient is therefore")
print("mostly a measure of how cleanly that reference can be scored, not of how much it matters - which is a warning against")
print("reading any single reference's estimate as biology, including the ones this project has quoted for YidC, NuoL and NuoN.")
print(f"\nBelow the floor the disattenuation is division by near-zero and the column is meaningless there: the apparent extremes")
print(f"({', '.join(ok.nlargest(1, 'disattenuated').ref)} at reliability {ok.nlargest(1, 'disattenuated').reliability.iloc[0]:.3f}, "
      f"{', '.join(ok.nsmallest(1, 'disattenuated').ref)} at {ok.nsmallest(1, 'disattenuated').reliability.iloc[0]:.3f}) are arithmetic, not evidence.")
print(f"\nRankable references (reliability >= {FLOOR:.2f}), by disattenuated coefficient:")
for y in rk.itertuples():
    print(f"   {y.ref} {y.name[:44]:44s} ({y.pfam}) rel {y.reliability:.3f}  raw {y.beta_per_sd:+.4f}  disatt {y.disattenuated:+.4f}  n={y.n}")
print(f"   correlation with the PC1 loading among these: {rk[['disattenuated', 'pc1_loading']].corr().iloc[0, 1]:+.3f}")
# These two correlations are the most-quoted figures from this script and were printed but never stored, so the
# synthesis quoted them from a terminal rather than from a file. scripts/69 caught that; they are now rows.
_summary = pd.DataFrame([
    {"ref": "SUMMARY", "name": "raw per-reference coefficient against that reference's reliability",
     "n": len(ok), "beta_per_sd": float(ok[["beta_per_sd", "reliability"]].corr().iloc[0, 1])},
    {"ref": "SUMMARY", "name": f"disattenuated coefficient against PC1 loading, reliability >= {FLOOR:.2f}",
     "n": len(rk), "beta_per_sd": float(rk[["disattenuated", "pc1_loading"]].corr().iloc[0, 1])},
])
pd.concat([R, _summary], ignore_index=True).to_csv(D / "reference_ranking.csv", index=False)
print(f"   (both correlations appended to {D / 'reference_ranking.csv'} as SUMMARY rows)")
d = R[R.pfam == DESIGN_PFAM]
if len(d):
    x = d.iloc[0]
    print(f"\nThe design's target family {DESIGN_PFAM} is represented in this panel by {x['ref']} ({x['name']}), and it is one of the")
    print(f"weakest references in it: reliability {x.reliability:.3f} - its exposed-periplasmic measure barely tracks the rest of the")
    print(f"panel at all - and a raw coefficient of {x.beta_per_sd:+.4f}, which is zero. The same Pfam family scored from NuoL and NuoN as")
    print("reference structures gave +0.076 and +0.088 per point in the original panel. The difference is the REFERENCE STRUCTURE, not")
    print("the family, so this is a reference-structure sensitivity of the kind the two-SecY-species check (p 0.42) was too narrow to")
    print("detect, and it bears directly on the design: the experiment's predicted effect was scaled from NuoL and NuoN specifically.")
