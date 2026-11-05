"""When the same Pfam family is scored from different reference structures, do the measurements agree? -> DATA_DIR/reference_structure_sensitivity.csv
Script 54 found that PF00361 - the NuoL/NuoM/NuoN family the charge-reversal experiment targets - carries the signal when scored from NuoL or NuoN as the reference structure
and nothing at all when scored from HyfB. That is a sensitivity to a modelling choice nobody has treated as a choice: the structural class of every position comes from ONE
AlphaFold model of ONE Escherichia coli protein, and every genome's homolog is then scored by alignment to it. The only check ever run on this was two SecY reference species
giving +0.042 and +0.071 with a difference of p 0.42, which compared two close homologues of the same protein and could not detect the problem.
This generalises the check. Three families in the project's two panels have more than one reference structure:
  PF00361  NuoL (P33607), NuoN (P0AFF0), HyfB (P23482)        - the design's target family, and the case that raised the question
  PF00344  SecY E. coli (P0AGA2), SecY B. subtilis (P16336)   - the pair already reported, re-run here on the same footing
  PF00115  CyoB appears in both panels under the same accession, so it is a consistency check on the pipeline rather than on structures
For each family: the correlation between the per-genome measurements taken from each reference structure (the same protein in the same genome, scored against different
reference topologies, so a low correlation is a measurement problem and nothing else), each reference's own coefficient against pH_min on the SAME genomes with composition
and source controlled, and whether the references differ by more than their standard errors allow.
A high correlation would mean the structural annotation transfers and the choice of reference is immaterial. A low one means the project's per-protein numbers are specific to
the reference structures chosen, and the aggregate result is carried by whatever survives averaging over that noise.
Usage: python3 scripts/55_reference_structure_sensitivity.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from itertools import combinations
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
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
# label -> (which results file holds it, pfam, human name)
GROUPS = {
    "PF00361 (NuoL/NuoM/NuoN, the design's target)": [("NuoL", "orig", "NuoL, P33607"), ("NuoN", "orig", "NuoN, P0AFF0"), ("P23482", "wide", "HyfB, P23482")],
    "PF00344 (SecY)": [("SecY", "orig", "SecY E. coli, P0AGA2"), ("SecY_Bsub", "orig", "SecY B. subtilis, P16336"), ("P0AGA2", "wide", "SecY E. coli again, via the wide panel")],
    "PF00115 (cytochrome bo3 subunit 1)": [("CyoB", "orig", "CyoB, P0ABI8"), ("P0ABI8", "wide", "CyoB again, via the wide panel")],
}
orig = pd.read_csv(D / "structure_panel.csv")
wide = pd.read_csv(D / "structure_panel_wide.csv")
reflen_o = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size() if (ROOT / "results/structure_panel_classes.csv").exists() else None
reflen_w = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
for src, tab in (("orig", orig), ("wide", wide)):
    tab.drop_duplicates("ncbi_assembly", inplace=True)
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))


def measure(label, src):
    """Per-genome exposed-periplasmic acidic excess from one reference structure, aligned to B's row order."""
    tab, reflen = (orig, reflen_o) if src == "orig" else (wide, reflen_w)
    c = f"{label}|peri_exp|n"
    if c not in tab:
        return None
    m = tab.set_index("ncbi_assembly")
    need = m[c].reindex(B.ncbi_assembly)
    floor = 8.0
    if reflen is not None and (label, "peri_exp") in reflen:
        floor = max(8, 0.5 * reflen[(label, "peri_exp")])
    num = (m[f"{label}|peri_exp|D"] + m[f"{label}|peri_exp|E"] - m[f"{label}|peri_exp|K"] - m[f"{label}|peri_exp|R"]).reindex(B.ncbi_assembly)
    x = (num / need).where(need >= floor)
    return pd.Series(x.to_numpy(), index=B.index)


rows = []
for fam, members in GROUPS.items():
    M = {}
    for label, src, nice in members:
        x = measure(label, src)
        if x is not None and x.notna().sum() > 200:
            M[nice] = x
    if len(M) < 2:
        print(f"{fam}: fewer than two usable reference structures, skipped\n")
        continue
    X = pd.DataFrame(M)
    keep = (X.notna().all(axis=1) & B.ph_min.notna()).to_numpy()
    n = int(keep.sum())
    print(f"{fam}: {len(M)} reference structures, {n} genomes scorable by all of them")
    print("   correlations between the measurements (same protein, same genome, different reference topology)")
    Xc = X[keep]
    for a, b_ in combinations(M, 2):
        r = float(Xc[a].corr(Xc[b_]))
        rows.append({"family": fam, "test": "correlation", "a": a, "b": b_, "n": n, "value": r})
        print(f"      {a[:34]:34s} vs {b_[:34]:34s}  r = {r:+.3f}")
    g = B[keep]
    C = Vcv[np.ix_(keep, keep)]
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    Xd = np.column_stack([np.ones(n)] + [g[c].to_numpy(float) for c in use])
    y = g.ph_min.to_numpy(float)
    lam = ph.pgls(y, Xd, C)["lambda"]
    Cn = C / np.mean(np.diag(C))
    L = np.linalg.cholesky(Cn * lam + np.diag(np.diag(Cn)) * (1 - lam) + 1e-10 * np.eye(n))
    Xs = solve_triangular(L, Xd, lower=True)
    Q, _ = np.linalg.qr(Xs)
    ys = solve_triangular(L, y, lower=True)
    yt = ys - Q @ (Q.T @ ys)
    dfree = n - Xs.shape[1] - 1
    print(f"   each reference's own coefficient on these same {n} genomes (per SD, composition and source controlled, lambda {lam:.3f})")
    est = {}
    for nice in M:
        v = Xc[nice].to_numpy(float)
        Fs = solve_triangular(L, (v - v.mean()) / v.std(), lower=True)
        Ft = Fs - Q @ (Q.T @ Fs)
        xx = float(Ft @ Ft)
        b = float(Ft @ yt) / xx
        s2 = (float(yt @ yt) - b * b * xx) / dfree
        se = np.sqrt(s2 / xx)
        est[nice] = (b, se)
        rows.append({"family": fam, "test": "coefficient", "a": nice, "n": n, "value": b, "se": se, "p": 2 * stats.t.sf(abs(b / se), dfree)})
        print(f"      {nice[:40]:40s} {b:+.4f} +/- {se:.4f} (p {2 * stats.t.sf(abs(b / se), dfree):.4f})")
    for a, b_ in combinations(est, 2):
        ba, sa = est[a]; bb, sb = est[b_]
        z = (ba - bb) / np.sqrt(sa ** 2 + sb ** 2)
        rows.append({"family": fam, "test": "difference", "a": a, "b": b_, "n": n, "value": ba - bb, "se": float(np.sqrt(sa ** 2 + sb ** 2)), "p": 2 * stats.norm.sf(abs(z))})
        print(f"      difference {a[:22]:22s} - {b_[:22]:22s} {ba - bb:+.4f} (z {z:+.2f}, p {2 * stats.norm.sf(abs(z)):.4f})"
              f"{'   <- they disagree' if 2 * stats.norm.sf(abs(z)) < 0.05 else ''}")
    print()

R = pd.DataFrame(rows)
R.to_csv(D / "reference_structure_sensitivity.csv", index=False)
print("wrote", D / "reference_structure_sensitivity.csv")
cr = R[R.test == "correlation"]
print(f"\nAcross all families, correlations between reference structures of the same Pfam run {cr.value.min():+.3f} to {cr.value.max():+.3f} "
      f"(median {cr.value.median():+.3f}).")
dd = R[(R.test == "difference") & (R.p < 0.05)]
print(f"{len(dd)} of {int((R.test == 'difference').sum())} pairs of reference structures give coefficients that differ by more than chance.")
