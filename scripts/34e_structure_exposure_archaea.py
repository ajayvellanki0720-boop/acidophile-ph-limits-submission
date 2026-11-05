"""Exposure test (34a/34b/34d) in Archaea (231 genomes on the archaeal marker tree). SecY (Sec61-alpha homolog) and EF-Tu (EF-1alpha) only: archaeal genomes carry no YidC-family hit aligning the -> DATA_DIR/structure_exposure_archaea_halophile_flag.csv
periplasmic domain. 'Periplasmic' here means the extracellular (pseudo-periplasmic, S-layer-facing) side of the reference topology. Models: total effect (all scored), total effect and linear / cubic
optimum on the optimum-matched subset (identical genomes). Halophiles = class Halobacteria (73 genomes; mean proteome pI 4.5-5.2 vs >= 5.74 in every other archaeon, so taxonomy and proteome signature agree); flagged in
structure_exposure_archaea_halophile_flag.csv. Usage: python3 scripts/34e_structure_exposure_archaea.py [--no-halophiles | --halophile-covariate | --halophiles-only] [--include-low-completeness]   (env DATA_DIR=results/archaea_ext PHYLO_DIR=results/phylo_arch_ext for the enlarged sample)"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from Bio import Phylo
from acidophile import phylo as ph
warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
CL = {"SecY": ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"], "EF-Tu": ["surface", "buried"]}
df = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "structure_exposure.csv"), on="ncbi_assembly")
tree = Phylo.read(T / "Archaea.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
A = df[(df.domain == "Archaea") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
A["halophile"] = A["class"].eq("Halobacteria")
if "completeness_lt90" in A and "--include-low-completeness" not in sys.argv:
    A = A[~A.completeness_lt90.fillna(False).astype(bool)].reset_index(drop=True)  # CheckM < 90 % kept only in the sensitivity run
LOWC = "--include-low-completeness" in sys.argv
NOHALO = "--no-halophiles" in sys.argv
HALOCOV = "--halophile-covariate" in sys.argv
HALONLY = "--halophiles-only" in sys.argv
if HALONLY:
    A = A[A.halophile].reset_index(drop=True)
    print(len(A), "halophilic archaea only")
if NOHALO:
    A[["ncbi_assembly", "genus", "class", "halophile"]].to_csv(D / "structure_exposure_archaea_halophile_flag.csv", index=False)
    A = A[~A.halophile].reset_index(drop=True)
    print(len(A), "non-halophile archaea")
SRC = ("thermobase", "acidb", "literature_auto", "bacdive")
for o in SRC: A[f"src_{o}"] = (A.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(A.ncbi_assembly))
z = (A.ph_opt - A.ph_opt.mean()) / A.ph_opt.std(); A["o1"], A["o2"], A["o3"] = z, z ** 2, z ** 3
rows = []
for ref, classes in CL.items():
    for c in classes:
        n = A[f"{ref}|{c}|n"]; A[f"{ref}_{c}"] = ((A[f"{ref}|{c}|D"] + A[f"{ref}|{c}|E"] - A[f"{ref}|{c}|K"] - A[f"{ref}|{c}|R"]) / n).where(n >= 15) * 100
    for c in classes:
        col = f"{ref}_{c}"
        for lab, opt, need_opt in (("total, all scored", [], False), ("total, optimum-matched", [], True), ("linear optimum", ["o1"], True), ("cubic optimum", ["o1", "o2", "o3"], True)):
            m = A[col].notna().to_numpy() & (A.ph_opt.notna().to_numpy() if need_opt else True)
            g = A[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
            src = [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
            if HALOCOV and g.halophile.nunique() > 1:
                g["halo"] = g.halophile.astype(float); src = src + ["halo"]
            X = np.column_stack([np.ones(len(g))] + [g[k] for k in opt + src])
            b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, g[[col]].to_numpy(float))
            rows.append({"ref": ref, "class": c, "model": lab, "n": len(g), "sd_pts": g[col].std(), "beta_per_pt": b[0], "se": se[0], "p": p[0]})
name = ("no_halophiles" if NOHALO else "halophiles_only" if HALONLY else "halo_covariate" if HALOCOV else "")
R = pd.DataFrame(rows); R.to_csv(D / (("lowc_" if LOWC else "") + "structure_exposure_archaea" + (f"_{name}" if name else "") + ".csv"), index=False)
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))
