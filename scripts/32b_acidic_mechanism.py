"""Where does the acidic-residue association live? PGLS (total effect: source indicators + feature(s)) on the broadened bacterial set, with the charge composition computed -> DATA_DIR/acidic_mechanism.csv
separately for membrane, secreted-like (signal-peptide heuristic), cytoplasmic, conserved-core (29 ribosomal-protein markers: the SAME protein families in every genome) and
non-core soluble proteins. Predictions: surface exposure -> membrane/secreted carry it; a proteome-wide shift -> the conserved core carries it too.
Usage: python3 scripts/32b_acidic_mechanism.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id")
      .merge(pd.read_csv(D / "compartment_charge.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["log_genome"] = np.log10(B.genome_size_bp)
C = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float)
z = lambda v: (np.asarray(v, float) - np.nanmean(v)) / np.nanstd(v)
CLS = ["all", "membrane", "secreted", "cytoplasmic", "core", "noncore_soluble"]
print(f"bacteria {len(B)} | share of residues: " + ", ".join(f"{c} {B[f'{c}_share_residues'].mean():.2f}" for c in CLS[1:4]) + f" | core (ribosomal) proteins per genome median {B.core_n_proteins.median():.0f}")
print("mean acidic excess (D+E-K-R per residue, %): " + ", ".join(f"{c} {100 * B[f'{c}_acidic_excess'].mean():+.2f}" for c in CLS))
print("corr of acidic excess between classes (core vs cytoplasmic %.2f, core vs membrane %.2f, membrane vs secreted %.2f)" % (B.core_acidic_excess.corr(B.cytoplasmic_acidic_excess), B.core_acidic_excess.corr(B.membrane_acidic_excess), B.membrane_acidic_excess.corr(B.secreted_acidic_excess)))


def fit(cols, extra=(), mask=None, with_opt=False):
    m = np.ones(len(B), bool) if mask is None else np.asarray(mask)
    g = B[m].reset_index(drop=True); Cm = C[np.ix_(m, m)]
    cov = (["ph_opt"] if with_opt else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1] + list(extra)
    X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in cov] + [z(g[c]) for c in cols])
    r = ph.pgls(g.ph_min.to_numpy(float), X, Cm); k = X.shape[1] - len(cols)
    return {c: (r["beta"][k + i], r["se"][k + i], r["p"][k + i]) for i, c in enumerate(cols)}, len(g), r["lambda"]


rows = []
print("\nA. each protein class alone (beta = pH_min units per SD of that class's acidic excess)")
for c in CLS:
    for lab, ex, wo, mk in (("total effect", (), False, None), ("+ GC + genome size", ("gc", "log_genome"), False, None), ("optimum as covariate", (), True, B.ph_opt.notna().to_numpy())):
        res, n, lam = fit([f"{c}_acidic_excess"], ex, mk, wo)
        b, se, p = res[f"{c}_acidic_excess"]; rows.append({"table": "A_single_class", "class": c, "model": lab, "n": n, "beta_per_sd": b, "se": se, "p": p})
A = pd.DataFrame(rows); print(A.pivot(index="class", columns="model", values="beta_per_sd").loc[CLS].round(3).to_string()); print(A.pivot(index="class", columns="model", values="p").loc[CLS].applymap(lambda x: f"{x:.0e}").to_string())

print("\nB. joint models")
for lab, cols in (("membrane + secreted + cytoplasmic", ["membrane_acidic_excess", "secreted_acidic_excess", "cytoplasmic_acidic_excess"]),
                  ("conserved core + non-core soluble", ["core_acidic_excess", "noncore_soluble_acidic_excess"]),
                  ("conserved core + membrane + secreted", ["core_acidic_excess", "membrane_acidic_excess", "secreted_acidic_excess"])):
    res, n, lam = fit(cols)
    print(f"  {lab}  (n={n}):  " + " | ".join(f"{c.replace('_acidic_excess', '')} {b:+.3f} ± {se:.3f} (p {p:.0e})" for c, (b, se, p) in res.items()))
    for c, (b, se, p) in res.items():
        rows.append({"table": "B_joint:" + lab, "class": c, "model": "joint", "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\nC. acidic (D+E) vs basic (K+R) content within the conserved core and within all proteins")
for c in ("core", "all", "membrane", "secreted", "cytoplasmic"):
    res, n, lam = fit([f"{c}_DE", f"{c}_KR"])
    print(f"  {c:12s} D+E {res[f'{c}_DE'][0]:+.3f} ± {res[f'{c}_DE'][1]:.3f} (p {res[f'{c}_DE'][2]:.0e}) | K+R {res[f'{c}_KR'][0]:+.3f} ± {res[f'{c}_KR'][1]:.3f} (p {res[f'{c}_KR'][2]:.0e})")
    for k, (b, se, p) in res.items():
        rows.append({"table": "C_DE_vs_KR", "class": k, "model": "joint DE+KR", "n": n, "beta_per_sd": b, "se": se, "p": p})

print("\nD. the conserved core by phylum (same proteins in every genome): pH per +1 percentage point of core acidic excess")
rows2 = []
for grp, gi in B.groupby("phylum"):
    if len(gi) < 40:
        continue
    m = (B.phylum == grp).to_numpy(); res, n, lam = fit(["core_acidic_excess"], mask=m); b, se, p = res["core_acidic_excess"]; sd = 100 * B[m].core_acidic_excess.std()
    rows2.append({"phylum": grp, "n": n, "core_SD_pts": sd, "pH_per_pt": b / sd, "se": se / sd, "p": p})
P = pd.DataFrame(rows2); print(P.round(3).to_string(index=False))
w = 1 / P.se ** 2; fe = (w * P.pH_per_pt).sum() / w.sum(); Q = (w * (P.pH_per_pt - fe) ** 2).sum(); k = len(P) - 1; tau2 = max(0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (P.se ** 2 + tau2); mu = (wr * P.pH_per_pt).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
print(f"  pooled {mu:+.3f} ± {sem:.3f} pH per percentage point (p {2 * stats.norm.sf(abs(mu / sem)):.1e}); I² {max(0, (Q - k) / Q):.0%}; positive in {(P.pH_per_pt > 0).sum()}/{len(P)} phyla")
pd.concat([pd.DataFrame(rows), P.assign(table="D_core_by_phylum")], ignore_index=True).to_csv(D / "acidic_mechanism.csv", index=False)
