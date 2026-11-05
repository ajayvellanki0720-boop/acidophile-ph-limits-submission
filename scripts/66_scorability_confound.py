"""Is the central result partly a measure of how well a genome aligns to an E. coli reference? -> DATA_DIR/scorability_confound.csv
Every per-genome value in this project is produced by aligning that genome's protein to one of 32 *E. coli* reference structures and reading residues at
the reference's exposed periplasmic positions. A genome is scored on a reference only if enough of those positions align, so how many of the 32 a genome
can be scored on, and how completely, is a property of how divergent it is from *E. coli*. The panel's coverage has been discussed throughout as a
property of the REFERENCES - which ones are prevalent enough to include - and never as a property of the GENOMES, which is where it would be a confound.
It is a confound of exactly the kind this project has been wrong about before, and the raw numbers say it is not hypothetical: across the bacteria the
number of references passing their floor correlates -0.136 with pH_min and mean coverage fraction correlates -0.187. Less alignable genomes have lower
reported minima. Acidophiles are phylogenetically clustered, so divergence from *E. coli* and acid tolerance are expected to travel together.
PGLS does not dispose of this. The lambda-scaled Brownian covariance absorbs variance shared by relatives, but distance from one particular reference
organism is a specific direction on the tree rather than autocorrelation among the tips, and nothing in the model knows about it.
Five parts.
  1. THE MEASURES, and what they are associated with: the count of references passing their floor, mean coverage fraction, and - since the composite
     needs 3 references - the composite's own missingness.
  2. DOES SCORABILITY PREDICT THE OUTCOME on its own, under the same PGLS the headline uses.
  3. IS IT CORRELATED WITH THE PREDICTOR, which is the other leg a confound needs. A scorability effect on pH_min matters only to the extent that
     scorability also moves the composite.
  4. THE HEADLINE WITH IT CONTROLLED, both measures, together and separately, and restricted to the well-scored genomes where it barely varies.
  5. THE SPECIFICITY CHECK, which is the part that can actually exonerate the claim. If scorability is driving the result it is an artefact of the
     alignment pipeline and should hit the matched cytoplasmic and transmembrane classes the same way, because they are scored by the same alignments
     against the same references. A nuisance that moves all classes equally cannot produce a periplasm-specific coefficient.
Usage: python3 scripts/66_scorability_confound.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
SRCC = [f"src_{o}" for o in SRC]
rows = []


def class_composite(cls):
    """The project's composite for one position class, plus the per-genome scorability of that class."""
    cols, ncols, covcols = {}, {}, {}
    for r in MEM:
        c = f"{r}|{cls}|n"
        if c not in wide.columns or (r, cls) not in reflen:
            continue
        n_ = wide[c].reindex(B.ncbi_assembly)
        floor = max(8, 0.5 * reflen[(r, cls)])
        num = (wide[f"{r}|{cls}|D"] + wide[f"{r}|{cls}|E"] - wide[f"{r}|{cls}|K"] - wide[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly)
        x = pd.Series((num / n_).where(n_ >= floor).to_numpy(), index=B.index)
        cols[r] = (x - x.mean()) / x.std()
        ncols[r] = pd.Series((n_ >= floor).to_numpy(), index=B.index)
        covcols[r] = pd.Series((n_ / reflen[(r, cls)]).to_numpy(), index=B.index)
    Z = pd.DataFrame(cols)
    return (Z.mean(axis=1).where(Z.notna().sum(axis=1) >= 3),
            pd.DataFrame(ncols).sum(axis=1).astype(float),
            pd.DataFrame(covcols).mean(axis=1))


pr, scor, cov = class_composite("peri_exp")
B["scor"] = scor
B["cov"] = cov
print("1. the measures")
print(f"   references passing their floor, per genome: mean {scor.mean():.2f}, SD {scor.std():.2f}, range {int(scor.min())}-{int(scor.max())}")
print(f"   mean coverage fraction: mean {cov.mean():.3f}, SD {cov.std():.3f}")
print(f"   genomes with no composite at all (fewer than 3 references): {int(pr.isna().sum())} of {len(B)}")
print(f"   raw correlation with pH_min: scorability {B.scor.corr(B.ph_min):+.3f}, coverage {B['cov'].corr(B.ph_min):+.3f}")
rows.append({"part": "measures", "mean_scorable": float(scor.mean()), "sd_scorable": float(scor.std()),
             "mean_coverage": float(cov.mean()), "r_scor_phmin": float(B.scor.corr(B.ph_min)),
             "r_cov_phmin": float(B["cov"].corr(B.ph_min)), "n_no_composite": int(pr.isna().sum())})


def fit(x, covs, mask=None):
    m = (x.notna() & B.ph_min.notna() & B[covs].notna().all(axis=1) if covs else x.notna() & B.ph_min.notna())
    if mask is not None:
        m = m & mask
    m = m.to_numpy()
    g = B[m]
    use = [c for c in covs + SRCC if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    b, se, p, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, x[m].to_numpy(float)[:, None])
    return float(b[0]), float(se[0]), float(p[0]), int(m.sum())


print("\n2. does scorability predict pH_min on its own, under the headline's model")
for name, v in (("references passing the floor", B.scor), ("mean coverage fraction", B["cov"])):
    b, se, p, n = fit(v, ["gw"])
    print(f"   {name:30s} n={n:4d}  {b:+.4f} +/- {se:.4f} per SD (p {p:.2e})")
    rows.append({"part": "scorability as predictor", "measure": name, "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\n3. is it correlated with the composite, which is the other leg a confound needs")
ok = pr.notna()
print(f"   composite against scorability r {pr[ok].corr(B.scor[ok]):+.3f};  against coverage r {pr[ok].corr(B['cov'][ok]):+.3f}")
rows.append({"part": "predictor correlation", "r_composite_scor": float(pr[ok].corr(B.scor[ok])),
             "r_composite_cov": float(pr[ok].corr(B["cov"][ok]))})

print("\n4. the headline with scorability controlled")
SPECS = [("composition only (the published model)", ["gw"], None),
         ("+ references passing the floor", ["gw", "scor"], None),
         ("+ mean coverage fraction", ["gw", "cov"], None),
         ("+ both", ["gw", "scor", "cov"], None),
         ("well-scored genomes only (>= 28 of 32)", ["gw"], B.scor >= 28),
         ("well-scored genomes, + both", ["gw", "scor", "cov"], B.scor >= 28)]
for name, covs, mask in SPECS:
    b, se, p, n = fit(pr, covs, mask)
    print(f"   {name:40s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {p:.2e})")
    rows.append({"part": "headline", "model": name, "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\n5. the specificity check: scorability is a property of the alignment, so it should hit every class alike")
print(f"   {'class':12s} {'published':>20s} {'+ scorability + coverage':>26s}  {'change':>8s}")
for cls in ("peri_exp", "cyto_exp", "TM", "peri_bur", "cyto_bur"):
    x, sc, cv = class_composite(cls)
    if x is None or x.notna().sum() < 200:
        print(f"   {cls:12s} not computable")
        continue
    B["scor_c"], B["cov_c"] = sc, cv
    b0, se0, p0, n0 = fit(x, ["gw"])
    b1, se1, p1, n1 = fit(x, ["gw", "scor_c", "cov_c"])
    print(f"   {cls:12s} {b0:+9.4f} (p {p0:7.1e}) {b1:+13.4f} (p {p1:7.1e})  {b1 - b0:+8.4f}")
    rows.append({"part": "specificity", "class": cls, "n": n0, "beta_published": b0, "p_published": p0,
                 "beta_controlled": b1, "p_controlled": p1, "change": b1 - b0})

R = pd.DataFrame(rows)
R.to_csv(D / "scorability_confound.csv", index=False)
print("\nwrote", D / "scorability_confound.csv")
print("\nReading. Scorability is a real correlate of the outcome, so it belongs in the record whatever it does to the coefficient. What")
print("decides whether it threatens the claim is part 5: a pipeline artefact scores every position class through the same alignments, so")
print("it cannot produce a coefficient in the periplasmic class and not in the matched cytoplasmic one. If controlling it leaves the")
print("periplasmic estimate standing and the controls null, the confound is disposed of rather than merely survived.")
