"""The archaeal test, redone on the wide panel -> results/archaea_ext2/archaea_wide_panel.csv
The earlier archaeal test (34e) could use only SecY as a membrane reference, because the hand-picked panel's other membrane proteins have no archaeal homolog that aligns the
periplasmic domain, and it was underpowered to the point of being uninformative. Of the 32 rule-selected references in the wide panel, 21 are found in at least 60 % of these
archaea, which is a different proposition: the composite can be built from many proteins rather than one.
References are kept when at least MINCOV of archaea can be scored for the class, so the membrane and soluble composites rest on comparable coverage. 'Periplasmic' means the
extracellular side of the reference topology; archaea have no periplasm, so for them the class is the S-layer-facing surface.
Halophiles (class Halobacteria) are handled three ways, as in 34e: excluded, included with an indicator, and alone - their proteomes are extremely acidic and their pH minima
are confounded with alkaliphily, so they cannot simply be pooled in.
Usage: python3 scripts/43_archaea_wide_panel.py"""
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
D = ROOT / "results/archaea_ext2"
T = ROOT / "results/phylo_arch_ext2"
MINCOV = 0.60
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
panel = pd.concat([pd.read_csv(ROOT / "results/broad2/structure_panel_wide.csv"), pd.read_csv(D / "structure_panel_wide_new.csv")], ignore_index=True)
pheno = pd.read_csv(D / "phenotypes_all.csv")
df = pheno.merge(panel, on="ncbi_assembly")
tree = Phylo.read(T / "Archaea.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
A = df[df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
A["halophile"] = A["class"].eq("Halobacteria")
if "completeness_lt90" in A:
    A = A[~A.completeness_lt90.fillna(False).astype(bool)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "literature_auto", "bacdive")
for o in SRC:
    A[f"src_{o}"] = (A.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(A.ncbi_assembly))
z = lambda s: (s - s.mean()) / s.std()
A["gw"] = A.acidic_excess * 100
print(f"{len(A)} archaea on the tree ({int(A.halophile.sum())} halophiles)\n")


def excess(ref, cls):
    c = f"{ref}|{cls}|n"
    if c not in A or (ref, cls) not in reflen:
        return None
    n = A[c]
    return ((A[f"{ref}|{cls}|D"] + A[f"{ref}|{cls}|E"] - A[f"{ref}|{cls}|K"] - A[f"{ref}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(ref, cls)]))


def composite(names, cls):
    cols, used = [], []
    for r in names:
        x = excess(r, cls)
        if x is None or x.notna().mean() < MINCOV:
            continue
        cols.append(z(x)); used.append(r)
    if not cols:
        return None, []
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= max(3, len(cols) // 2)), used


rows = []
for cls, names, label in (("peri_exp", MEM, "exposed extracellular (membrane)"), ("cyto_exp", MEM, "exposed cytoplasmic (membrane)"),
                          ("peri_bur", MEM, "buried extracellular (membrane)"), ("TM", MEM, "transmembrane"), ("surface", SOL, "soluble surface")):
    pr, used = composite(names, cls)
    if pr is None:
        continue
    print(f"{label}: {len(used)} references pass {MINCOV:.0%} archaeal coverage")
    for sub, mask, cov in (("halophiles excluded", ~A.halophile, ["gw"]),
                           ("all, halophile indicator", pd.Series(True, index=A.index), ["gw", "halo"]),
                           ("halophiles only", A.halophile, ["gw"])):
        m = (pr.notna() & mask).to_numpy()
        if m.sum() < 60:
            continue
        g = A[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
        g["halo"] = g.halophile.astype(float)
        use = [c for c in cov + [f"src_{o}" for o in SRC] if c in g and g[c].nunique() > 1]
        X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
        b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, pr[m].to_numpy(float)[:, None])
        rows.append({"class": label, "refs": len(used), "subset": sub, "n": len(g), "beta_per_sd": b[0], "se": se[0], "p": p[0]})
        print(f"   {sub:26s} n={len(g):3d}  {b[0]:+.4f} +/- {se[0]:.4f} (p {p[0]:.3f})")
R = pd.DataFrame(rows); R.to_csv(D / "archaea_wide_panel.csv", index=False)
print("\nBacterial comparison for the same class and model: +0.185 per SD (composition controlled).")
e = R[(R["class"].str.startswith("exposed extracellular")) & (R.subset == "halophiles excluded")]
if len(e):
    b, se = float(e.beta_per_sd.iloc[0]), float(e.se.iloc[0])
    print(f"Archaeal estimate, halophiles excluded: {b:+.4f} +/- {se:.4f}; "
          f"difference from bacteria {b - 0.185:+.4f} ({abs(b - 0.185) / se:.1f} archaeal SE); "
          f"detectable effect at 80% power {2.8 * se:.3f}")
    print(f"one-sided p for the bacterial value being too large: {stats.norm.sf((0.185 - b) / se):.3f}")
