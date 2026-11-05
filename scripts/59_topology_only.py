"""Does the result need AlphaFold, or only UniProt topology? -> DATA_DIR/topology_only.csv
Script 55 showed that per-genome measurements barely transfer between reference structures of the same Pfam family (median correlation +0.126), and the structural half of the
annotation is the part that depends on a predicted model: which positions face the periplasm comes from UniProt's curated topological domains, but whether a position is
solvent-exposed comes from solvent-accessible surface area computed on a single-chain AlphaFold model, filtered by per-residue pLDDT. If the claim can be restated using
topology alone it does not depend on any structure prediction, and the reference-structure noise stops being a threat to it. If it cannot, then AlphaFold is load-bearing and
that noise is a real worry.
The test is a decomposition that costs nothing new, because the per-class tallies are already computed. For each reference, the periplasmic positions split into exposed and
buried, and the two can be recombined:
  periplasmic, ALL positions     - topology only, no structure prediction used beyond deciding which side a position is on
  periplasmic, exposed only      - the published claim
  periplasmic, buried only       - the complement
and the same three on the cytoplasmic side as the matched control. Recombination is done on the COUNTS, not by averaging the two classes' excesses, so a reference whose
periplasmic positions are mostly buried is weighted accordingly rather than counted as half exposed.
The informative comparisons are (a) whether periplasmic-all survives at all, which decides whether AlphaFold is needed; (b) how much of the exposed-only effect it retains,
which measures how much the exposure filter is worth; and (c) whether the cytoplasmic side stays null under the same treatment, since a topology-only result would be worth
little if it appeared on both sides.
Usage: python3 scripts/59_topology_only.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
cls_tab = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv")
reflen = cls_tab.groupby(["ref", "class"]).size()
wide = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
W = wide.set_index("ncbi_assembly")
idx = B.ncbi_assembly

# how the reference positions themselves divide, which sets the ceiling on what topology alone can carry
pos = cls_tab[cls_tab["class"].isin(["peri_exp", "peri_bur", "cyto_exp", "cyto_bur"])].groupby("class").size()
print("reference positions by class across the 32 membrane references:")
for c in ("peri_exp", "peri_bur", "cyto_exp", "cyto_bur"):
    print(f"   {c:9s} {int(pos.get(c, 0)):5d}")
pe, pb = int(pos.get("peri_exp", 0)), int(pos.get("peri_bur", 0))
print(f"   so {pe / (pe + pb):.0%} of periplasmic reference positions are exposed\n")


def tally(label, classes):
    """Acidic excess over the union of several classes, recombined on counts rather than on per-class excesses."""
    num = pd.Series(0.0, index=idx)
    den = pd.Series(0.0, index=idx)
    ok = pd.Series(False, index=idx)
    floor = 0.0
    for cls in classes:
        c = f"{label}|{cls}|n"
        if c not in W or (label, cls) not in reflen:
            continue
        n_ = W[c].reindex(idx).fillna(0.0)
        d = (W[f"{label}|{cls}|D"] + W[f"{label}|{cls}|E"] - W[f"{label}|{cls}|K"] - W[f"{label}|{cls}|R"]).reindex(idx).fillna(0.0)
        num = num + d
        den = den + n_
        ok = ok | W[c].reindex(idx).notna()
        floor += 0.5 * reflen[(label, cls)]
    if not ok.any():
        return None
    x = (num / den).where(ok & (den >= max(8, floor)))
    return pd.Series(x.to_numpy(), index=B.index)


def composite(classes):
    cols = {}
    for r in MEM:
        x = tally(r, classes)
        if x is None or x.notna().sum() < 200:
            continue
        cols[r] = (x - x.mean()) / x.std()
    Z = pd.DataFrame(cols)
    return Z.mean(axis=1).where(Z.notna().sum(axis=1) >= 3), len(cols)


def fit(x):
    m = (x.notna() & B.ph_min.notna()).to_numpy()
    g = B[m]
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    b, se, p, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, x[m].to_numpy(float)[:, None])
    return float(b[0]), float(se[0]), float(p[0]), int(m.sum())


SPECS = [
    ("periplasmic, ALL positions (topology only)", ["peri_exp", "peri_bur"]),
    ("periplasmic, exposed only (the published claim)", ["peri_exp"]),
    ("periplasmic, buried only", ["peri_bur"]),
    ("cytoplasmic, ALL positions (topology only)", ["cyto_exp", "cyto_bur"]),
    ("cytoplasmic, exposed only", ["cyto_exp"]),
    ("cytoplasmic, buried only", ["cyto_bur"]),
]
rows, est = [], {}
print("composite estimates, genome-wide composition controlled (per SD)")
for label, classes in SPECS:
    pr, k = composite(classes)
    if pr is None or pr.notna().sum() < 200:
        print(f"   {label:48s} not computable")
        continue
    b, se, p, n = fit(pr)
    est[label] = (b, se)
    rows.append({"predictor": label, "classes": "+".join(classes), "references": k, "n": n, "beta_per_sd": b, "se": se, "p": p})
    print(f"   {label:48s} refs {k:2d}  n={n:4d}  {b:+.4f} +/- {se:.4f} (p {p:.2e})")

R = pd.DataFrame(rows)
R.to_csv(D / "topology_only.csv", index=False)
print("\nwrote", D / "topology_only.csv")
topo = est.get("periplasmic, ALL positions (topology only)")
expo = est.get("periplasmic, exposed only (the published claim)")
ctopo = est.get("cytoplasmic, ALL positions (topology only)")
if topo and expo:
    print(f"\nTopology alone retains {topo[0] / expo[0]:.0%} of the exposed-only estimate ({topo[0]:+.4f} against {expo[0]:+.4f}).")
    z = (expo[0] - topo[0]) / np.sqrt(expo[1] ** 2 + topo[1] ** 2)
    print(f"The exposure filter adds {expo[0] - topo[0]:+.4f}, which is z {z:+.2f} - not an independent test, since the two composites share their exposed positions.")
    if ctopo:
        print(f"The matched cytoplasmic side under the same treatment: {ctopo[0]:+.4f} +/- {ctopo[1]:.4f}, so the topology-only result "
              f"{'is specific to the periplasmic side' if abs(ctopo[0]) < abs(topo[0]) / 2 else 'appears on both sides and is NOT specific'}.")
    print("\nReading: if the topology-only row stands on its own, the claim does not depend on AlphaFold and the reference-structure")
    print("noise in scripts/55 threatens the exposure refinement rather than the result.")
