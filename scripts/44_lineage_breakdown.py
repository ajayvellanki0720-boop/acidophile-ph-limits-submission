"""Does the wide-panel result hold lineage by lineage? -> DATA_DIR/lineage_breakdown.csv
The histidine lead died on exactly this test: it was carried by two phyla and absent in two others. Four things here, all on the 32-protein exposed-periplasmic composite with
genome-wide composition controlled:
  1. each phylum with at least 60 species on its own, then random-effects pooling and the heterogeneity across them;
  2. leave-one-phylum-out, to check no single lineage carries the pooled estimate;
  3. within genera (genus-mean-centred, permuted inside genus), a finer scale than the family test - if the association only appears between genera this is where it fails;
  4. the within-family estimate expressed per standard deviation of the RAW composite as well as of the centred one, because the synthesis table mixes the two scales and they
     are not the same number (the centred composite has the smaller spread).
Usage: python3 scripts/44_lineage_breakdown.py [n_permutations]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
z = lambda s: (s - s.mean()) / s.std()
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC]


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


def fit(mask, pr):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    use = [c for c in ["gw"] + SRCC if g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, pr[mask].to_numpy(float)[:, None])
    return b[0], se[0], p[0], len(g)


def pool(bs, ses):
    bs, ses = np.asarray(bs), np.asarray(ses)
    w = 1 / ses ** 2; fe = (w * bs).sum() / w.sum(); Q = (w * (bs - fe) ** 2).sum(); k = len(bs) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (ses ** 2 + tau2); mu = (wr * bs).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    i2 = max(0.0, (Q - k) / Q) if k > 0 and Q > 0 else 0.0
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), i2, (1 - stats.chi2.cdf(Q, k)) if k > 0 else np.nan


rows = []
pr = PR["exposed periplasmic"]
big = [p for p, n in B[pr.notna()].phylum.value_counts().items() if n >= 60]
print(f"{len(B)} bacteria; phyla with >= 60 scored species: {', '.join(big)}\n")
print("1. each phylum on its own (beta per SD, composition controlled)")
bs, ses = [], []
for p_ in big:
    m = (pr.notna() & (B.phylum == p_)).to_numpy()
    b, se, pv, n = fit(m, pr)
    bs.append(b); ses.append(se)
    rows.append({"test": "per phylum", "group": p_, "n": n, "beta_per_sd": b, "se": se, "p": pv})
    print(f"   {p_:22s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {pv:.4f})")
mu, sem, pv, i2, ph_q = pool(bs, ses)
print(f"   pooled {mu:+.4f} +/- {sem:.4f} (p {pv:.1e}); positive in {sum(1 for x in bs if x > 0)}/{len(bs)}; I2 {i2:.0%}; heterogeneity p {ph_q:.2f}")
rows.append({"test": "per phylum pooled", "group": "random effects", "n": int(pr.notna().sum()), "beta_per_sd": mu, "se": sem, "p": pv})

print("\n2. leave-one-phylum-out")
for p_ in big:
    m = (pr.notna() & (B.phylum != p_)).to_numpy()
    b, se, pv, n = fit(m, pr)
    rows.append({"test": "leave-one-phylum-out", "group": f"without {p_}", "n": n, "beta_per_sd": b, "se": se, "p": pv})
    print(f"   without {p_:22s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {pv:.2e})")

print("\n3. within genera and within families, permuted inside the group")
rng = np.random.default_rng(20261026)
for name, pr_ in PR.items():
    for level, minn in (("genus", 3), ("family", 4)):
        d = B.assign(x=pr_)[[level, "ph_min", "x", "gw"]].dropna()
        d = d[d[level].isin(d[level].value_counts()[lambda s: s >= minn].index)]
        if len(d) < 200:
            continue
        A = np.column_stack([np.ones(len(d)), d.gw])
        xr = d.x - A @ np.linalg.lstsq(A, d.x, rcond=None)[0]
        yr = d.ph_min - A @ np.linalg.lstsq(A, d.ph_min, rcond=None)[0]
        dd = d.assign(x=xr, y=yr)
        c = dd.groupby(level)[["y", "x"]].transform(lambda s: s - s.mean())
        obs = np.polyfit(c.x, c.y, 1)[0]
        gidx = [np.flatnonzero((dd[level] == f).to_numpy()) for f in dd[level].unique()]
        xv, yv = c.x.to_numpy(), c.y.to_numpy(); null = np.empty(NPERM)
        for i in range(NPERM):
            xp = xv.copy()
            for ix in gidx:
                xp[ix] = rng.permutation(xp[ix])
            null[i] = np.polyfit(xp, yv, 1)[0]
        pv = (1 + np.sum(np.abs(null) >= abs(obs))) / (NPERM + 1)
        sd_raw, sd_cen = pr_.std(), c.x.std()
        rows.append({"test": f"within {level}", "group": name, "n": len(dd), "beta_per_sd": obs * sd_cen, "se": null.std() * sd_cen, "p": pv,
                     "beta_per_sd_raw_scale": obs * sd_raw, "groups": dd[level].nunique()})
        print(f"   {name:20s} within {level:7s} n={len(dd):4d} in {dd[level].nunique():3d} groups: "
              f"{obs * sd_cen:+.4f} per SD of centred composite, {obs * sd_raw:+.4f} on the raw-composite scale (permutation p {pv:.4f})")
pd.DataFrame(rows).to_csv(D / "lineage_breakdown.csv", index=False)
print("\nwrote", D / "lineage_breakdown.csv")
