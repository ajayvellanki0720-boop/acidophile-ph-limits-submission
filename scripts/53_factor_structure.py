"""Is there one latent surface-charge quantity, or is the composite just an average? -> DATA_DIR/factor_structure.csv
Script 52 left a sharp question. The whole interpretation of the composite - and the explanation of why 32 references beat 7 - rests on the references measuring one common
quantity plus independent noise. But their mean pairwise correlation is only +0.09 to +0.11, which is weak enough that the one-factor model deserves to be checked directly
rather than inferred from the fact that averaging helps.
Three tests on the per-reference exposed-periplasmic measures:
  1. HOW MANY FACTORS. Eigenvalues of the correlation matrix against a parallel-analysis null: each column is permuted independently, which destroys the covariance between
     references while preserving each one's own distribution, and the test is how many observed eigenvalues exceed the 95th percentile of the permuted ones. One factor clearly
     above the null and the rest inside it is what the interpretation needs.
  2. DO THE LOADINGS AGREE IN SIGN. If the references measure a common quantity, the first component should load them all the same way. Mixed signs would mean the unweighted
     mean is cancelling real structure rather than averaging noise, and that the composite works for some other reason.
  3. DOES THE FACTOR BEAT THE AVERAGE. The first component is the maximum-variance weighting of the references; if one latent quantity drives the association it should predict
     the lower limit at least as well as the unweighted mean. The comparison is fair because the component is estimated WITHOUT the outcome, so there is no leakage - but the
     loadings are estimated on the complete-case subset and then applied to every genome, which is where the two predictors are compared.
Usage: python3 scripts/53_factor_structure.py [n_permutations]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 500
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100


def zmat(names, cls):
    out = {}
    for r in names:
        c = f"{r}|{cls}|n"
        if c not in B or (r, cls) not in reflen:
            continue
        n_ = B[c]
        x = ((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n_).where(n_ >= max(8, 0.5 * reflen[(r, cls)]))
        out[r] = (x - x.mean()) / x.std()
    return pd.DataFrame(out)


Z = zmat(MEM, "peri_exp")
cov = Z.notna().mean()
names = [r for r in Z.columns if cov[r] >= 0.80]
cc = Z[names].notna().all(axis=1).to_numpy()
M = Z.loc[cc, names].to_numpy(float)
k = len(names)
print(f"{k} references at >=80% coverage, {cc.sum()} complete-case genomes")
R = np.corrcoef(M, rowvar=False)
tri = R[np.triu_indices(k, 1)]
print(f"pairwise correlations: mean {tri.mean():+.3f}, range {tri.min():+.3f} to {tri.max():+.3f}, {int((tri > 0).sum())}/{len(tri)} positive\n")

print("1. how many factors, against a parallel-analysis null")
ev = np.linalg.eigvalsh(R)[::-1]
rng = np.random.default_rng(53)
pev = np.empty((NPERM, k))
for i in range(NPERM):
    P = np.column_stack([rng.permutation(M[:, j]) for j in range(k)])
    pev[i] = np.linalg.eigvalsh(np.corrcoef(P, rowvar=False))[::-1]
thr = np.percentile(pev, 95, axis=0)
nfac = int(np.sum(ev > thr))
rows = []
for j in range(min(6, k)):
    rows.append({"test": "eigenvalue", "component": j + 1, "value": float(ev[j]), "null_p95": float(thr[j]),
                 "variance_explained": float(ev[j] / k)})
    print(f"   component {j + 1}: eigenvalue {ev[j]:.3f} against null 95th percentile {thr[j]:.3f}   "
          f"({ev[j] / k:.1%} of variance){'  <- above null' if ev[j] > thr[j] else ''}")
print(f"   {nfac} component(s) exceed the null; the first explains {ev[0] / k:.1%} of the variance in {k} references")

print("\n2. first-component loadings")
w = np.linalg.eigh(R)[1][:, -1]
if w.sum() < 0:
    w = -w
same = int((w > 0).sum())
equal = 1 / np.sqrt(k)  # the loading every reference would carry if they contributed equally
print(f"   {same}/{k} loadings positive; range {w.min():+.3f} to {w.max():+.3f}; equal weighting would give {equal:+.3f} each, "
      f"and the loadings' spread is {w.std() / abs(w.mean()):.2f} times their mean")
for r_, wi in sorted(zip(names, w), key=lambda z: -z[1])[:3]:
    print(f"      highest: {r_} {wi:+.3f}")
for r_, wi in sorted(zip(names, w), key=lambda z: z[1])[:3]:
    print(f"      lowest:  {r_} {wi:+.3f}")
rows += [{"test": "loading", "component": 1, "ref": r_, "value": float(wi)} for r_, wi in zip(names, w)]

print("\n3. does the factor predict the limit better than the unweighted mean?")
# loadings from the complete-case subset, applied to every genome; both predictors use whatever references a genome has
Zn = Z[names]
mean_pr = Zn.mean(axis=1).where(Zn.notna().sum(axis=1) >= 3)
Wd = pd.DataFrame(np.tile(w, (len(Zn), 1)), columns=names, index=Zn.index).where(Zn.notna())
pc_pr = (Zn * Wd).sum(axis=1, min_count=3) / Wd.abs().sum(axis=1)
pc_pr = pc_pr.where(Zn.notna().sum(axis=1) >= 3)
soluble = zmat(SOL, "surface")
sol_pr = soluble.mean(axis=1).where(soluble.notna().sum(axis=1) >= 3)
keep = (mean_pr.notna() & pc_pr.notna() & B.ph_min.notna()).to_numpy()
g = B[keep].reset_index(drop=True)
C = ph.vcv_fast(tree, list(B.ncbi_assembly))[np.ix_(keep, keep)]
use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in use])
y = g.ph_min.to_numpy(float)
F = np.column_stack([mean_pr[keep].to_numpy(float), pc_pr[keep].to_numpy(float)])
b, se, p, lam = ph.pgls_scan(C, y, X, F)
print(f"   n={len(g)}, lambda {lam:.3f}")
for lbl, i in (("unweighted mean of the references", 0), ("first principal component", 1)):
    rows.append({"test": "prediction", "group": lbl, "n": len(g), "value": float(b[i]), "se": float(se[i]), "p": float(p[i])})
    print(f"   {lbl:34s} {b[i]:+.4f} +/- {se[i]:.4f} per SD (p {p[i]:.2e})")
print(f"   correlation between the two predictors: {np.corrcoef(F[:, 0], F[:, 1])[0, 1]:+.4f}")
# both together, to see whether the component carries anything the mean does not
fit = ph.pgls(y, np.column_stack([X, (F - F.mean(0)) / F.std(0)]), C, lam=lam)
print(f"   both in one model: mean {fit['beta'][-2]:+.4f} (p {fit['p'][-2]:.3f}), component {fit['beta'][-1]:+.4f} (p {fit['p'][-1]:.3f})")
rows.append({"test": "prediction", "group": "both, mean coefficient", "n": len(g), "value": float(fit["beta"][-2]), "se": float(fit["se"][-2]), "p": float(fit["p"][-2])})
rows.append({"test": "prediction", "group": "both, component coefficient", "n": len(g), "value": float(fit["beta"][-1]), "se": float(fit["se"][-1]), "p": float(fit["p"][-1])})
pd.DataFrame(rows).to_csv(D / "factor_structure.csv", index=False)
print("\nwrote", D / "factor_structure.csv")
print(f"\nReading: {nfac} factor(s) above the parallel-analysis null, the first explaining {ev[0] / k:.0%} of the variance, "
      f"{same}/{k} loadings of one sign, and the component {'does' if abs(b[1]) > abs(b[0]) else 'does not'} beat the unweighted mean.")
