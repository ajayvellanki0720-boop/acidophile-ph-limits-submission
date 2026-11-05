"""Does the result lean on genomes that are not from the strain whose pH was measured? -> DATA_DIR/genome_provenance.csv
The phenotype table has carried two provenance columns since it was built and neither has ever been tested, because both looked like bookkeeping.
`assembly_type_material` records whether the genome sequenced is from the species' type strain, and `assembly_level` how finished the assembly is. The
first turns out to matter a great deal more than bookkeeping, for a reason that is specific to this project's shape.
Only 36 of 2,771 bacterial genomes are not from type material - but their mean pH_min is 2.62 against 5.80 for the rest, and they make up 29.3% of the
acid-tolerant tail at pH_min 3 or below (24 of 82). They are the canonical acidophiles: Acidithiobacillus, Leptospirillum, Acidiphilium,
Methylacidiphilum, Sulfobacillus. The reason is mundane - extreme acidophiles are hard to culture, so their sequenced genome is more often a convenient
strain rather than the type strain - and the consequence is not: the association is one-sided and lives almost entirely in that tail (scripts/51, 63),
so 36 rows with a possible genome-to-phenotype mismatch sit exactly where the result comes from. A composite computed on one strain against a pH
measured on another is measurement error in the predictor with no reason to be independent of the outcome.
Four parts.
  1. THE PROVENANCE VARIABLES and their association with pH_min, fitted under the headline's own model rather than eyeballed from means.
  2. THE CRITICAL REFIT: drop the 36 and see what is left. Their leverage is far out of proportion to their number, so a large drop here would be the
     most serious finding of the week; a small one disposes of the worry.
  3. ASSEMBLY LEVEL, which is the milder of the two. A fragmented assembly loses proteins, which should show up as scorability (tested in scripts/66)
     rather than as anything new, so this is mostly a check that it does.
  4. SPECIFICITY, as always: a provenance artefact is scored through the same alignments for every position class and cannot be periplasm-specific.
A thing this test cannot do. It asks whether the result depends on those 36 rows, not whether their phenotypes are right. Verifying the match between
each sequenced strain and the strain whose growth was measured needs the two strain designations compared by hand against the species description, and
for species whose description is paywalled that is not available here.
Usage: python3 scripts/67_genome_provenance.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
wide = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly").set_index("ncbi_assembly")
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
B["nontype"] = (~B.assembly_type_material.astype(bool)).astype(float)
for lv in ("Contig", "Scaffold", "Chromosome"):          # Complete Genome is the reference level
    B[f"lvl_{lv}"] = (B.assembly_level == lv).astype(float)
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
SRCC = [f"src_{o}" for o in SRC]
rows = []


def class_composite(cls):
    cols = {}
    for r in MEM:
        c = f"{r}|{cls}|n"
        if c not in wide.columns or (r, cls) not in reflen:
            continue
        n_ = wide[c].reindex(B.ncbi_assembly)
        num = (wide[f"{r}|{cls}|D"] + wide[f"{r}|{cls}|E"] - wide[f"{r}|{cls}|K"] - wide[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly)
        x = pd.Series((num / n_).where(n_ >= max(8, 0.5 * reflen[(r, cls)])).to_numpy(), index=B.index)
        cols[r] = (x - x.mean()) / x.std()
    Z = pd.DataFrame(cols)
    return Z.mean(axis=1).where(Z.notna().sum(axis=1) >= 3)


def fit(x, covs, mask=None):
    m = x.notna() & B.ph_min.notna() & B[covs].notna().all(axis=1) if covs else x.notna() & B.ph_min.notna()
    if mask is not None:
        m = m & mask
    m = m.to_numpy()
    g = B[m]
    use = [c for c in covs + SRCC if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    b, se, p, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, x[m].to_numpy(float)[:, None])
    return float(b[0]), float(se[0]), float(p[0]), int(m.sum())


pr = class_composite("peri_exp")
nt = B.nontype.astype(bool)
print("1. the provenance variables")
print(f"   non-type-material genomes: {int(nt.sum())} of {len(B)}; mean pH_min {B.ph_min[nt].mean():.2f} against {B.ph_min[~nt].mean():.2f}")
tail = B.ph_min <= 3
print(f"   they are {nt[tail].mean():.1%} of the acid-tolerant tail (pH_min <= 3), which is {int(nt[tail].sum())} of {int(tail.sum())} rows")
print(f"   their composite: mean {pr[nt].mean():+.3f} against {pr[~nt].mean():+.3f} for the rest")
print(f"   assembly level: " + ", ".join(f"{k} {v}" for k, v in B.assembly_level.value_counts().items()))
rows.append({"part": "description", "n_nontype": int(nt.sum()), "n_total": len(B), "mean_phmin_nontype": float(B.ph_min[nt].mean()),
             "mean_phmin_type": float(B.ph_min[~nt].mean()), "tail_share_nontype": float(nt[tail].mean()),
             "composite_nontype": float(pr[nt].mean()), "composite_type": float(pr[~nt].mean())})
for name, v in (("non-type indicator", B.nontype), ("Contig assembly", B.lvl_Contig), ("Scaffold assembly", B.lvl_Scaffold)):
    if v.nunique() < 2:
        continue
    b, se, p, n = fit(v, ["gw"])
    print(f"   {name:22s} as a predictor of pH_min: {b:+.4f} +/- {se:.4f} per SD (p {p:.2e})")
    rows.append({"part": "provenance as predictor", "measure": name, "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\n2. the critical refit: the headline without the 36")
SPECS = [("the published model", ["gw"], None),
         ("+ non-type indicator", ["gw", "nontype"], None),
         ("DROP the 36 non-type genomes", ["gw"], ~nt),
         ("drop the 36, acid-tolerant tail only (pH_min <= 4)", ["gw"], (~nt) & (B.ph_min <= 4)),
         ("the published model, tail only (pH_min <= 4)", ["gw"], B.ph_min <= 4)]
for name, covs, mask in SPECS:
    b, se, p, n = fit(pr, covs, mask)
    print(f"   {name:50s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {p:.2e})")
    rows.append({"part": "refit", "model": name, "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\n3. assembly level, which should reduce to scorability")
lvls = [c for c in ("lvl_Contig", "lvl_Scaffold", "lvl_Chromosome") if B[c].nunique() > 1]
for name, covs in (("the published model", ["gw"]), ("+ assembly level indicators", ["gw"] + lvls),
                   ("+ assembly level and non-type", ["gw", "nontype"] + lvls)):
    b, se, p, n = fit(pr, covs)
    print(f"   {name:34s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {p:.2e})")
    rows.append({"part": "assembly level", "model": name, "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\n4. specificity, with provenance controlled")
print(f"   {'class':12s} {'published':>20s} {'+ provenance':>20s} {'drop the 36':>20s}")
for cls in ("peri_exp", "cyto_exp", "TM"):
    x = class_composite(cls)
    if x is None or x.notna().sum() < 200:
        continue
    b0, s0, p0, n0 = fit(x, ["gw"])
    b1, s1, p1, n1 = fit(x, ["gw", "nontype"] + lvls)
    b2, s2, p2, n2 = fit(x, ["gw"], ~nt)
    print(f"   {cls:12s} {b0:+9.4f} (p {p0:6.1e}) {b1:+9.4f} (p {p1:6.1e}) {b2:+9.4f} (p {p2:6.1e})")
    rows.append({"part": "specificity", "class": cls, "beta_published": b0, "p_published": p0,
                 "beta_provenance": b1, "p_provenance": p1, "beta_dropped": b2, "p_dropped": p2})

R = pd.DataFrame(rows)
R.to_csv(D / "genome_provenance.csv", index=False)
print("\nwrote", D / "genome_provenance.csv")
print("\nReading. 36 rows at 1.3% of the sample but 29% of the tail that carries the association is the kind of leverage that decides a")
print("result, so the number to read is the refit without them. If it holds, the worry is disposed of and the provenance columns become")
print("a recorded caveat rather than a threat. If it does not, the finding was always resting on genomes that may not belong to the")
print("strains whose growth was measured, and no amount of phylogenetic machinery would have shown it.")
