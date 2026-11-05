"""Envelope type as a moderator across ALL sampled phyla -> DATA_DIR/envelope_metaregression.csv
Script 47 compared two monoderm phyla against two diderm ones, so envelope type was confounded with phylogeny. Breaking that confound inside Bacillota is not possible here: the
dataset holds only 7 Negativicutes, the diderm class of that phylum. The alternative is to use more phyla. Ten have at least 12 scored species - three monoderm (Bacillota,
Actinomycetota, Chloroflexota) and seven diderm (Pseudomonadota, Bacteroidota, Thermodesulfobacteriota, Campylobacterota, Acidobacteriota, Verrucomicrobiota, Planctomycetota) -
which supports a random-effects meta-regression of each phylum's own estimate on its envelope type rather than a two-against-two contrast.
Deinococcota and Thermotogota are held out of the moderator test and reported separately: both have atypical envelopes (Deinococcota stain Gram-positive but carry an outer
membrane, Thermotogota have an outer sheath), so assigning them either way would be a judgement call that drives the result.
The small phyla give noisy estimates - a within-phylum PGLS on 16 to 34 species is weak - and precision weighting handles that, but the moderator test still rests mostly on the
four large phyla, so it is a strengthening of 47 rather than an independent confirmation.
Usage: python3 scripts/48_envelope_metaregression.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph
warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MONODERM = {"Bacillota", "Actinomycetota", "Chloroflexota"}
DIDERM = {"Pseudomonadota", "Bacteroidota", "Thermodesulfobacteriota", "Campylobacterota", "Acidobacteriota", "Verrucomicrobiota", "Planctomycetota"}
ATYPICAL = {"Deinococcota", "Thermotogota", "Aquificota"}
MINN = 12
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
z = lambda s: (s - s.mean()) / s.std()
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in ("thermobase", "acidb", "bacdive", "literature_auto")]
def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen: continue
        n = B[f"{r}|{cls}|n"]
        cols.append(z(((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
def fit(mask, pr):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    use = [c for c in ["gw"] + SRCC if g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, pr[mask].to_numpy(float)[:, None])
    return b[0], se[0], p[0], len(g)
PRED = {"exposed periplasmic": composite(MEM, "peri_exp"), "exposed cytoplasmic": composite(MEM, "cyto_exp"), "soluble surface": composite(SOL, "surface")}
rows = []
for name, pr in PRED.items():
    counts = B[pr.notna()].phylum.value_counts()
    use = [p_ for p_, n in counts.items() if n >= MINN and p_ in (MONODERM | DIDERM | ATYPICAL)]
    est = []
    print(f"=== {name}")
    for p_ in use:
        m = (pr.notna() & (B.phylum == p_)).to_numpy()
        if m.sum() < MINN: continue
        b, se, pv, n = fit(m, pr)
        env = "monoderm" if p_ in MONODERM else ("diderm" if p_ in DIDERM else "atypical")
        est.append((p_, env, b, se, n))
        rows.append({"predictor": name, "phylum": p_, "envelope": env, "n": n, "beta_per_sd": b, "se": se, "p": pv})
        print(f"   {p_:26s} {env:9s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {pv:.3f})")
    # random-effects meta-regression on envelope, atypical phyla excluded
    e = [(b, se, 1.0 if env == "monoderm" else 0.0) for _, env, b, se, _ in est if env in ("monoderm", "diderm")]
    if len(e) >= 4:
        bs = np.array([x[0] for x in e]); ses = np.array([x[1] for x in e]); mono = np.array([x[2] for x in e])
        w = 1 / ses ** 2
        X = np.column_stack([np.ones(len(bs)), mono])
        beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], bs * np.sqrt(w), rcond=None)[0]
        resid = bs - X @ beta
        Q = float((w * resid ** 2).sum()); kdf = len(bs) - 2
        tau2 = max(0.0, (Q - kdf) / (w.sum() - np.trace(np.linalg.pinv(X.T @ (X * w[:, None])) @ (X.T @ (X * (w ** 2)[:, None]))))) if kdf > 0 else 0.0
        wr = 1 / (ses ** 2 + tau2)
        V = np.linalg.pinv(X.T @ (X * wr[:, None]))
        beta = V @ (X.T @ (wr * bs))
        se_b = np.sqrt(np.diag(V))
        pv = stats.norm.sf(beta[1] / se_b[1])
        print(f"   meta-regression over {len(bs)} phyla ({int(mono.sum())} monoderm): diderm intercept {beta[0]:+.4f} +/- {se_b[0]:.4f}, "
              f"monoderm increment {beta[1]:+.4f} +/- {se_b[1]:.4f}, one-sided p {pv:.4f}; tau2 {tau2:.4f}")
        rows.append({"predictor": name, "phylum": f"meta-regression ({len(bs)} phyla)", "envelope": "monoderm increment", "n": int(sum(x[4] for x in est)),
                     "beta_per_sd": beta[1], "se": se_b[1], "p": pv})
    print()
pd.DataFrame(rows).to_csv(D / "envelope_metaregression.csv", index=False)
