"""Is p = 1e-6 real, or is it what a phylogenetically autocorrelated predictor produces by chance? -> DATA_DIR/null_calibration.csv
Every p-value in this project comes from PGLS on a neighbour-joining tree with Pagel's lambda fitted by maximum likelihood. Two things could make those p-values too small, and
neither has been tested: the predictor is itself strongly phylogenetically autocorrelated (closely related genomes have similar surface composition), and the tree is wrong in
ways a bootstrap of its own alignment cannot see. If a smooth-on-the-tree predictor is tested against a response whose covariance the model gets wrong, the test can reject far
more often than 5 % of the time with no effect present at all.
This is a parametric bootstrap. The response is SIMULATED with no effect of the predictor - drawn from the null PGLS model's own multivariate normal, mean X*beta_null and
covariance sigma2 * V(lambda_gen) - while the predictor keeps its real values and therefore its real phylogenetic structure. The test is then run exactly as the project runs it.
Three generating lambdas, to separate sampling noise from model misspecification:
  lambda_gen = lambda_hat  the model is correct; the rejection rate should be the nominal one, and anything above it is the test's own fault;
  lambda_gen = 1           pure Brownian motion, more phylogenetic signal than the fitted model allows - the realistic misspecification, since lambda is fitted on the null model
                           and the tree is a distance tree;
  lambda_gen = 0           no phylogenetic signal, the opposite error.
The test always uses lambda_hat, held fixed across simulations. Fixing it is a concession: refitting lambda per simulation needs a Cholesky factorisation of a 1,657 x 1,657
matrix about twenty times per draw, which is hours rather than minutes. Holding it fixed removes one source of variability from the null, so these rates are if anything
optimistic about the test, which is the safe direction for a check of this kind.
Reported for both the exposed-periplasmic composite and the soluble-surface control: rejection rates at 0.05 and 0.001, and where the REAL |t| falls in the simulated null.
Usage: python3 scripts/50_null_calibration.py [n_simulations]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
from scipy.linalg import solve_triangular
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
NSIM = int(sys.argv[1]) if len(sys.argv) > 1 else 500
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
z = lambda s: (s - s.mean()) / s.std()


def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen:
            continue
        n = B[f"{r}|{cls}|n"]
        cols.append(z(((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)


PR = {"exposed periplasmic": composite(MEM, "peri_exp"), "soluble surface": composite(SOL, "surface")}
keep = (PR["exposed periplasmic"].notna() & PR["soluble surface"].notna() & B.ph_min.notna()).to_numpy()
g = B[keep].reset_index(drop=True)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))[np.ix_(keep, keep)]
n = len(g)
use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
X = np.column_stack([np.ones(n)] + [g[c].to_numpy(float) for c in use])
y = g.ph_min.to_numpy(float)
XS = {k: v[keep].to_numpy(float) for k, v in PR.items()}
Cn = Call / np.mean(np.diag(Call))
print(f"{n} bacteria with both composites and a minimum; design {X.shape[1]} columns ({', '.join(use)})\n")


def chol(lam):
    V = Cn * lam + np.diag(np.diag(Cn)) * (1 - lam)
    return np.linalg.cholesky(V + 1e-10 * np.eye(n))


null = ph.pgls(y, X, Call)
lam_hat, s2_hat, beta_null = null["lambda"], null["sigma2"], null["beta"]
print(f"null model (covariates only): lambda_hat {lam_hat:.3f}, residual sigma2 {s2_hat:.4f}")

# One factorisation per lambda. The test side is built once at lambda_hat: whitened design, its orthonormal basis, and each predictor's residualised whitened column.
L_test = chol(lam_hat)
Xs = solve_triangular(L_test, X, lower=True)
Q, _ = np.linalg.qr(Xs)
dfree = n - Xs.shape[1] - 1
FT = {}
for name, x in XS.items():
    Fz = (x - x.mean()) / x.std()
    Fs = solve_triangular(L_test, Fz, lower=True)
    FT[name] = Fs - Q @ (Q.T @ Fs)


def tstat(yv):
    """|t| for each predictor, from pgls_scan's own algebra at lambda_hat, reusing the factorisation."""
    ys = solve_triangular(L_test, yv, lower=True)
    yt = ys - Q @ (Q.T @ ys)
    out = {}
    for name, Ft in FT.items():
        xx = float(Ft @ Ft)
        b = float(Ft @ yt) / xx
        s2 = (float(yt @ yt) - b * b * xx) / dfree
        out[name] = (b, abs(b) / np.sqrt(s2 / xx))
    return out


obs = tstat(y)
for name, (b, t) in obs.items():
    print(f"   real data: {name:20s} beta {b:+.4f} per SD, |t| {t:.2f}, p {2 * stats.t.sf(t, dfree):.2e}")
# agreement check against the function the project actually calls
chk, _, chkp, _ = ph.pgls_scan(Call, y, X, np.column_stack([XS[k] for k in XS]), lam=lam_hat)
print(f"   pgls_scan agreement: betas {chk[0]:+.4f}, {chk[1]:+.4f} (max difference {max(abs(chk[i] - list(obs.values())[i][0]) for i in range(2)):.2e})\n")

rows = []
for lam_gen in (lam_hat, 1.0, 0.0):
    Lg = chol(lam_gen)
    rng = np.random.default_rng(50)
    sim = {k: np.empty(NSIM) for k in XS}
    simb = {k: np.empty(NSIM) for k in XS}
    mu = X @ beta_null
    sd = np.sqrt(s2_hat)
    for i in range(NSIM):
        ysim = mu + sd * (Lg @ rng.standard_normal(n))
        for name, (b, t) in tstat(ysim).items():
            sim[name][i] = t
            simb[name][i] = b
    tag = "lambda_hat (model correct)" if lam_gen == lam_hat else f"lambda_gen = {lam_gen:.0f}"
    print(f"{NSIM} simulations with no effect, response generated at {tag}:")
    for name in XS:
        t = sim[name]
        r05 = float(np.mean(2 * stats.t.sf(t, dfree) < 0.05))
        r001 = float(np.mean(2 * stats.t.sf(t, dfree) < 0.001))
        exceed = float(np.mean(t >= obs[name][1]))
        rows.append({"lambda_gen": lam_gen, "predictor": name, "n_sim": NSIM, "reject_0.05": r05, "reject_0.001": r001,
                     "max_null_abs_t": float(t.max()), "real_abs_t": obs[name][1], "frac_null_above_real": exceed,
                     "null_beta_sd": float(simb[name].std()), "real_beta": obs[name][0]})
        print(f"   {name:20s} rejects at 0.05: {r05:5.1%} (nominal 5%)   at 0.001: {r001:5.1%} (nominal 0.1%)   "
              f"largest null |t| {t.max():5.2f} vs real {obs[name][1]:5.2f}; null betas SD {simb[name].std():.4f} vs real {obs[name][0]:+.4f}")
    print()
R = pd.DataFrame(rows)
R.to_csv(D / "null_calibration.csv", index=False)
print("wrote", D / "null_calibration.csv")
w = R[(R.lambda_gen == 1.0) & (R.predictor == "exposed periplasmic")].iloc[0]
infl = w["reject_0.05"] / 0.05
print(f"\nWorst case for the test (Brownian response, lambda fixed at the fitted value): the 5% test rejects {w['reject_0.05']:.1%} of the time, "
      f"an inflation of {infl:.1f}x. Even there the real |t| of {w['real_abs_t']:.1f} is {'beyond' if w['frac_null_above_real'] == 0 else 'inside'} "
      f"the whole simulated null (largest null |t| {w['max_null_abs_t']:.2f}).")
