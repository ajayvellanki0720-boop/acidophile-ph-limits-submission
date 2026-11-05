"""Redo the headline numbers with the factor score instead of the unweighted mean -> DATA_DIR/factor_score_headline.csv
Script 53 showed the first principal component of the per-reference measures predicts the lower limit better than the unweighted mean the project has used throughout
(+0.244 against +0.198 per standard deviation, and in a joint model the mean falls to -0.017 while the component holds +0.260). Every headline number in the synthesis was
computed with the mean, so each of them is a lower bound on what the better predictor gives, and the question is whether the difference matters anywhere it would change a
conclusion. Three places it might:
  1. The EFFECT SIZE with composition controlled and then with the growth optimum controlled - the two numbers quoted in the claim.
  2. The PRACTICAL VALUE, which is the conclusion most exposed to a weak predictor: the mean adds +0.0065 of out-of-sample R-squared over the growth optimum, and if that is
     mostly attenuation then the honest headline changes.
  3. The ACIDOPHILE RANKING, where the mean alone reaches AUC 0.785.
The loadings are estimated INSIDE each training fold, never on held-out genomes. Principal components do not use the outcome, so leakage is not the worry; the point is that
a predictor whose construction peeked at the test fold's covariance is not the predictor anyone could use prospectively, and the comparison should be of things that could
actually be deployed. The unweighted mean is recomputed in exactly the same folds so that the two differ only in the weighting.
Usage: python3 scripts/56_factor_score_headline.py [n_repeats]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
cols = {}
for r in MEM:
    c = f"{r}|peri_exp|n"
    if c not in B or (r, "peri_exp") not in reflen:
        continue
    n_ = B[c]
    x = ((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n_).where(n_ >= max(8, 0.5 * reflen[(r, "peri_exp")]))
    cols[r] = (x - x.mean()) / x.std()
Z = pd.DataFrame(cols)
names = [r for r in Z.columns if Z[r].notna().mean() >= 0.80]
print(f"{len(names)} references at >=80% coverage used for the weighting; {Z.shape[1]} scorable in all\n")


def pc_weights(Zsub):
    """First-component loadings from the complete cases of Zsub, signed so that they are mostly positive."""
    cc = Zsub.notna().all(axis=1)
    w = np.linalg.eigh(np.corrcoef(Zsub[cc].to_numpy(float), rowvar=False))[1][:, -1]
    return -w if w.sum() < 0 else w


def apply_weights(Zsub, w, minrefs=3):
    W = pd.DataFrame(np.tile(w, (len(Zsub), 1)), columns=Zsub.columns, index=Zsub.index).where(Zsub.notna())
    s = (Zsub * W).sum(axis=1, min_count=minrefs) / W.abs().sum(axis=1)
    return s.where(Zsub.notna().sum(axis=1) >= minrefs)


Zn = Z[names]
B["mean_pr"] = Zn.mean(axis=1).where(Zn.notna().sum(axis=1) >= 3)
B["pc_pr"] = apply_weights(Zn, pc_weights(Zn))
rows = []

print("1. effect size, same model as the synthesis quotes")
keep = (B.mean_pr.notna() & B.pc_pr.notna() & B.ph_min.notna()).to_numpy()
g = B[keep].reset_index(drop=True)
C = ph.vcv_fast(tree, list(B.ncbi_assembly))[np.ix_(keep, keep)]
y = g.ph_min.to_numpy(float)
SRCC = [f"src_{o}" for o in SRC]
for lbl, covs in (("composition controlled", ["gw"]), ("composition + growth optimum", ["gw", "ph_opt"])):
    m = (g[covs].notna().all(axis=1)).to_numpy()
    gg = g[m]
    use = [c for c in covs + SRCC if gg[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [gg[c].to_numpy(float) for c in use])
    F = np.column_stack([gg.mean_pr.to_numpy(float), gg.pc_pr.to_numpy(float)])
    b, se, p, lam = ph.pgls_scan(C[np.ix_(m, m)], y[m], X, F)
    rows.append({"test": "effect size", "model": lbl, "predictor": "unweighted mean", "n": int(m.sum()), "value": b[0], "se": se[0], "p": p[0]})
    rows.append({"test": "effect size", "model": lbl, "predictor": "factor score", "n": int(m.sum()), "value": b[1], "se": se[1], "p": p[1]})
    print(f"   {lbl:30s} n={int(m.sum()):4d}  mean {b[0]:+.4f} (p {p[0]:.1e})   factor {b[1]:+.4f} (p {p[1]:.1e})   ratio {b[1] / b[0]:.2f}x")

print("\n2. practical value: out-of-sample R-squared, 10 folds blocked by family, weighting learned inside each training fold")
cvrows = (B.mean_pr.notna() & B.ph_min.notna() & B.ph_opt.notna() & B.family.notna()).to_numpy()
d = B[cvrows].reset_index(drop=True)
Zd = Zn[cvrows].reset_index(drop=True)
print(f"   {len(d)} bacteria in {d.family.nunique()} families")
fams = d.family.unique()
rng = np.random.default_rng(56)
MODELS = {"optimum alone": None, "optimum + unweighted mean": "mean", "optimum + factor score": "pc"}
res = {k: [] for k in MODELS}
for rep in range(NREP):
    assign = {f: i % 10 for i, f in enumerate(rng.permutation(fams))}
    fold = d.family.map(assign).to_numpy()
    num = {k: 0.0 for k in MODELS}
    den = 0.0
    for k in range(10):
        tr, te = fold != k, fold == k
        if te.sum() < 5 or tr.sum() < 50:
            continue
        w = pc_weights(Zd[tr])                       # loadings from the training fold only
        pcs = apply_weights(Zd, w)
        preds = {"mean": d.mean_pr.to_numpy(float), "pc": pcs.to_numpy(float)}
        ytr, yte = d.loc[tr, "ph_min"].to_numpy(float), d.loc[te, "ph_min"].to_numpy(float)
        den += float(((yte - ytr.mean()) ** 2).sum())
        for name, key in MODELS.items():
            cs = [d.ph_opt.to_numpy(float)] + ([preds[key]] if key else [])
            A = np.column_stack(cs)
            mu, sg = A[tr].mean(0), A[tr].std(0); sg[sg == 0] = 1
            Atr = np.column_stack([np.ones(tr.sum()), (A[tr] - mu) / sg])
            Ate = np.column_stack([np.ones(te.sum()), (A[te] - mu) / sg])
            beta = np.linalg.lstsq(Atr, ytr, rcond=None)[0]
            num[name] += float(((yte - Ate @ beta) ** 2).sum())
    for name in MODELS:
        res[name].append(1 - num[name] / den)
for name in MODELS:
    v = np.array(res[name])
    rows.append({"test": "out-of-sample R2", "model": name, "n": len(d), "value": float(v.mean()), "se": float(v.std())})
    print(f"   {name:28s} R2 {v.mean():.4f} +/- {v.std():.4f}")
base = np.array(res["optimum alone"])
for name in ("optimum + unweighted mean", "optimum + factor score"):
    v = np.array(res[name])
    rows.append({"test": "R2 delta", "model": name, "n": len(d), "value": float((v - base).mean()), "se": float((v - base).std()),
                 "improved_folds": int((v > base).sum())})
    print(f"   delta from {name.split('+ ')[1]:18s} {(v - base).mean():+.4f}  (better in {int((v > base).sum())}/{len(v)})")

print("\n3. ranking extreme acidophiles (pH_min <= 3), same folds")
d["acid"] = (d.ph_min <= 3).astype(int)


def auc(lab, sc):
    pos, neg = sc[lab == 1], sc[lab == 0]
    r = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


for lbl, key in (("unweighted mean alone", "mean"), ("factor score alone", "pc")):
    a = []
    rng2 = np.random.default_rng(56)
    for rep in range(NREP):
        assign = {f: i % 10 for i, f in enumerate(rng2.permutation(fams))}
        fold = d.family.map(assign).to_numpy()
        pred = np.full(len(d), np.nan)
        for k in range(10):
            tr, te = fold != k, fold == k
            if te.sum() < 5:
                continue
            v = d.mean_pr.to_numpy(float) if key == "mean" else apply_weights(Zd, pc_weights(Zd[tr])).to_numpy(float)
            pred[te] = -v[te]  # a LOWER composite means more acid-tolerant, so negate to score acidophily
        m = ~np.isnan(pred)
        a.append(auc(d.acid.to_numpy()[m], pred[m]))
    a = np.array(a)
    rows.append({"test": "AUC", "model": lbl, "n": len(d), "value": float(a.mean()), "se": float(a.std())})
    print(f"   {lbl:24s} AUC {a.mean():.3f} +/- {a.std():.3f}")

# the two AUC GAINS are what the project actually quotes from this script, and they were printed only as a comparison
# in prose. correction 12's lesson is that a derived figure quoted anywhere has to be a stored row or it is unguarded.
_auc = {r["model"]: r["value"] for r in rows if r.get("test") == "AUC"}
_mean = next((v for k, v in _auc.items() if "unweighted" in k), None)
_fact = next((v for k, v in _auc.items() if "factor" in k), None)
if _mean is not None and _fact is not None:
    rows.append({"test": "AUC gain", "model": "factor score over the unweighted mean of the same 19 references",
                 "n": len(d), "value": float(_fact - _mean)})
    print(f"   gain from factor weighting: {_fact - _mean:+.3f} of AUC")
    try:
        _all32 = pd.read_csv(D / "shape_and_discrimination.csv")
        _b = _all32[(_all32.test == "AUC acidophile") & (_all32.group == "composite alone")].estimate
        if len(_b):
            rows.append({"test": "AUC gain", "model": "dropping the low-coverage references, 32 refs to the 19 high-coverage",
                         "n": len(d), "value": float(_mean - float(_b.iloc[0]))})
            print(f"   gain from dropping low-coverage references: {_mean - float(_b.iloc[0]):+.3f} of AUC "
                  f"(against the 32-reference mean's {float(_b.iloc[0]):.3f} from scripts/51)")
    except Exception as _e:
        print("   (32-reference AUC unavailable:", _e, ")")
pd.DataFrame(rows).to_csv(D / "factor_score_headline.csv", index=False)
print("\nwrote", D / "factor_score_headline.csv")
