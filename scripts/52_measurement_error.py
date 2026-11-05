"""Is the composite a noisy measure of one latent quantity, and how big is that quantity's effect? -> DATA_DIR/measurement_error.csv
Correction 4 in the synthesis asserts an explanation it never tested. The rule-selected panel of 32 references gives +0.185 per standard deviation against the hand-picked
panel of 7's +0.101, and I explained the difference by saying that averaging more references removes measurement error from the composite, so the standardised coefficient rises
towards the effect of whatever latent quantity the references all measure. That is a testable claim with a sharp prediction, in two parts.
  1. RELIABILITY. If the references measure one common quantity plus independent noise, splitting a panel in half gives two composites that correlate, and the Spearman-Brown
     formula turns that split-half correlation into the reliability of the full-length composite. Near 1 would mean the composite is already clean and there is nothing to
     correct; near 0 would mean it is mostly noise.
  2. SCALING. Classical measurement error attenuates a coefficient expressed per UNIT of the predictor by its reliability, but the coefficients here are per STANDARD DEVIATION
     of the predictor, and the standard deviation is itself inflated by the noise. Writing x = T + e with reliability rho = var(T)/var(x): the per-unit slope is attenuated to
     b*rho, and the per-SD slope is b*rho*sd(x) = b*sd(T)*sqrt(rho). So a per-SD coefficient disattenuates by dividing by sqrt(rho), NOT by rho. Dividing by rho over-corrects,
     by a further factor of 1/sqrt(rho). (This was got wrong on the first pass; `_check_disattenuation` below is the simulation that settled it, and it runs on every execution
     so the claim is never again taken on trust.) The prediction, then: beta_per_sd / sqrt(reliability) should be the SAME number at every panel size. If it keeps climbing with
     k, the references are not measuring one thing and the published explanation is wrong.
The design has to hold the sample fixed while k varies, or a change in which genomes are scorable gets read as a change in measurement. That is more constraining than it
looks: only ONE genome of 2,768 can be scored on all 32 references, because reference coverage ranges from 19% to 100%. So the test is run twice, trading panel size against
sample size, and the two have to agree for the conclusion to stand:
  A  the 19 references covering at least 80% of genomes, complete cases only (n about 440), k up to 16 - wide panel, small sample;
  B  the 10 references covering at least 85%, complete cases only (n about 1,440), k up to 8 - narrow panel, large sample.
Lambda is fitted once per configuration on the covariates-only model and held fixed, since it does not depend on which references the composite averages; that makes the fits
cheap enough to average over 25 random subsets at each k.
Usage: python3 scripts/52_measurement_error.py [n_subsets]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NSUB = int(sys.argv[1]) if len(sys.argv) > 1 else 25
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
cov = Z.notna().mean()
print(f"{Z.shape[1]} references scorable; coverage from {cov.min():.0%} to {cov.max():.0%}, only {int(Z.notna().all(axis=1).sum())} genome(s) scorable on all of them\n")

def _check_disattenuation(seed=0, n=400000, b=0.30):
    """A per-SD coefficient disattenuates by sqrt(reliability). Simulated here rather than asserted, because the first pass divided by reliability itself."""
    rng_ = np.random.default_rng(seed)
    out = []
    for rho in (0.25, 0.5, 0.8):
        T = rng_.normal(0, 1, n)
        x = T + rng_.normal(0, np.sqrt((1 - rho) / rho), n)
        y_ = b * T + rng_.normal(0, 1, n)
        bo = float(np.polyfit((x - x.mean()) / x.std(), y_, 1)[0])
        out.append((rho, bo, bo / rho, bo / np.sqrt(rho)))
    return out


print("sanity check on the correction itself (truth 0.3000):")
for rho, bo, by_rho, by_sqrt in _check_disattenuation():
    print(f"   reliability {rho:.2f}: observed {bo:.4f}  /reliability {by_rho:.4f} (wrong)  /sqrt(reliability) {by_sqrt:.4f} (right)")
print()

rows = []
for letter, thr, ks in (("A", 0.80, (1, 2, 4, 8, 12, 16)), ("B", 0.85, (1, 2, 4, 6, 8))):
    names = [r for r in Z.columns if cov[r] >= thr]
    tag = f"{letter}: {len(names)} refs at >={thr:.0%} coverage"
    keep = (Z[names].notna().all(axis=1) & B.ph_min.notna()).to_numpy()
    g = B[keep].reset_index(drop=True)
    Zc = Z.loc[keep, names].reset_index(drop=True)
    n = len(g)
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X = np.column_stack([np.ones(n)] + [g[c].to_numpy(float) for c in use])
    y = g.ph_min.to_numpy(float)
    C = Vcv[np.ix_(keep, keep)]
    lam = ph.pgls(y, X, C)["lambda"]
    Cn = C / np.mean(np.diag(C))
    L = np.linalg.cholesky(Cn * lam + np.diag(np.diag(Cn)) * (1 - lam) + 1e-10 * np.eye(n))
    Xs = solve_triangular(L, X, lower=True)
    Q, _ = np.linalg.qr(Xs)
    ys = solve_triangular(L, y, lower=True)
    yt = ys - Q @ (Q.T @ ys)
    dfree = n - Xs.shape[1] - 1
    tri = Zc.corr().to_numpy()[np.triu_indices(len(names), 1)]
    print(f"{tag}: {len(names)} references, {n} complete-case genomes, lambda {lam:.3f}, mean inter-reference r {np.nanmean(tri):+.3f}")

    def beta_per_sd(x):
        """pgls_scan's algebra at the fixed lambda, reusing one factorisation. Returns (beta per SD of x, t)."""
        xz = (x - x.mean()) / x.std()
        Fs = solve_triangular(L, xz, lower=True)
        Ft = Fs - Q @ (Q.T @ Fs)
        xx = float(Ft @ Ft)
        b = float(Ft @ yt) / xx
        s2 = (float(yt @ yt) - b * b * xx) / dfree
        return b, b / np.sqrt(s2 / xx)

    rng = np.random.default_rng(52)
    print(f"   {'k':>3}  {'beta/SD':>9}  {'spread':>7}  {'reliability':>11}  {'disattenuated':>13}")
    for k in ks:
        bs, rels = [], []
        nsub = 1 if k == len(names) else NSUB
        for _ in range(nsub):
            sub = list(rng.choice(names, size=k, replace=False)) if k < len(names) else list(names)
            b, t = beta_per_sd(Zc[sub].mean(axis=1).to_numpy(float))
            bs.append(b)
            if k >= 2:
                h = list(rng.permutation(sub))
                half = k // 2
                a, c2 = Zc[h[:half]].mean(axis=1), Zc[h[half:2 * half]].mean(axis=1)
                r = float(np.corrcoef(a, c2)[0, 1])
                # Spearman-Brown from two half-length forms up to the full length k
                rels.append(2 * r / (1 + r) if r > 0 else np.nan)
        bm, bsd = float(np.mean(bs)), float(np.std(bs))
        rel = float(np.nanmean(rels)) if rels else np.nan
        dis = bm / np.sqrt(rel) if rels and rel > 0 else np.nan  # per-SD coefficients disattenuate by sqrt(reliability); see _check_disattenuation
        rows.append({"config": tag, "k": k, "n_refs": len(names), "n": n, "n_subsets": nsub, "beta_per_sd": bm,
                     "beta_sd_across_subsets": bsd, "reliability": rel, "disattenuated": dis})
        print(f"   {k:>3}  {bm:>+9.4f}  {bsd:>7.4f}  {rel if rels else float('nan'):>11.3f}  {dis:>+13.4f}")
    full, tfull = beta_per_sd(Zc.mean(axis=1).to_numpy(float))
    print(f"   whole set of {len(names)}: {full:+.4f} per SD (|t| {abs(tfull):.2f}, p {2 * stats.t.sf(abs(tfull), dfree):.2e})\n")

R = pd.DataFrame(rows)
R.to_csv(D / "measurement_error.csv", index=False)
print("wrote", D / "measurement_error.csv")
print("\nPREDICTION: if the references measure one latent quantity plus independent noise, beta/sqrt(reliability) is the same at every k.")
for tag in R.config.unique():
    s = R[(R.config == tag) & R.disattenuated.notna()]
    if len(s) < 2:
        continue
    lo, hi = s.iloc[0], s.iloc[-1]
    print(f"  {tag}")
    print(f"     observed beta      {lo.beta_per_sd:+.4f} at k={int(lo.k)} -> {hi.beta_per_sd:+.4f} at k={int(hi.k)}  ({hi.beta_per_sd / lo.beta_per_sd:.2f}x)")
    print(f"     disattenuated      {lo.disattenuated:+.4f} at k={int(lo.k)} -> {hi.disattenuated:+.4f} at k={int(hi.k)}  ({hi.disattenuated / lo.disattenuated:.2f}x), "
          f"range over all k {s.disattenuated.min():+.4f} to {s.disattenuated.max():+.4f}")
