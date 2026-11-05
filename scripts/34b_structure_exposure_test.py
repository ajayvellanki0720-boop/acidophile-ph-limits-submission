"""Does the acidic-residue association with pH_min sit on exposed periplasmic surfaces? PGLS (bacteria) of pH_min on the acidic excess (D+E-K-R per residue) of each structural -> DATA_DIR/structure_exposure_test.csv
class of three conserved reference proteins (see 34a): membrane proteins YidC and SecY by topology/exposure, soluble EF-Tu (surface vs buried) as control.
Betas are per percentage point of acidic excess in that class so classes with different variances are comparable; the per-SD beta is also given. A class is scored only in
genomes that align >= MIN_N reference positions of it. Three models: total effect; optimum as covariate; total effect + genome-wide acidic excess (does the class add information?).
Usage: python3 scripts/34b_structure_exposure_test.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
MIN_N = 15
CLASSES = {"YidC": ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"], "SecY": ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"], "EF-Tu": ["surface", "buried"]}
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_exposure.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
rows = []
for ref, classes in CLASSES.items():
    for c in classes:
        n = B[f"{ref}|{c}|n"]
        x = (B[f"{ref}|{c}|D"] + B[f"{ref}|{c}|E"] - B[f"{ref}|{c}|K"] - B[f"{ref}|{c}|R"]) / n
        x = x.where(n >= MIN_N) * 100  # percentage points
        for lab in ("total effect", "optimum as covariate", "plus genome-wide acidic excess"):
            m = x.notna().to_numpy() & (B.ph_opt.notna().to_numpy() if lab == "optimum as covariate" else True)
            if m.sum() < 200:
                continue
            g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
            cov = (["ph_opt"] if lab == "optimum as covariate" else []) + (["acidic_excess"] if lab.startswith("plus") else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
            if lab.startswith("plus"):
                g["acidic_excess"] = g.acidic_excess * 100
            X = np.column_stack([np.ones(len(g))] + [g[k] for k in cov])
            xv = x[m].to_numpy(float)
            b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, xv[:, None])
            # BOTH SCALE COLUMNS BELOW ARE WRONG, found 2026-10-31 and kept with this note because the files are
            # already written and this script is superseded by the 36a panel. pgls_scan returns per SD of xv, so:
            #   `beta_per_pt` is per SD, not per percentage point (correction 15, divide by sd_pts to convert);
            #   `beta_per_sd` MULTIPLIES by the SD, which is correction 6's error verbatim - the result is neither
            #   per SD nor per unit but per unit times the SD squared, and it is quoted nowhere, which is how it
            #   survived here while correction 6 was being written about the same operation in other scripts.
            rows.append({"ref": ref, "class": c, "model": lab, "n": len(g), "mean_positions": float(n[m].mean()),
                         "sd_pts": xv.std(), "beta_per_pt": b[0], "se": se[0], "p": p[0],
                         "beta_per_sd": b[0] * xv.std(),
                         "scale_note": "beta_per_pt is per SD of the share; beta_per_sd is meaningless - see comment"})
R = pd.DataFrame(rows); R.to_csv(D / "structure_exposure_test.csv", index=False)
pd.set_option("display.width", 200)
for lab in ("total effect", "optimum as covariate", "plus genome-wide acidic excess"):
    print(f"\n{lab}: pH_min vs acidic excess of the structural class")
    print(R[R.model == lab][["ref", "class", "n", "mean_positions", "sd_pts", "beta_per_pt", "se", "beta_per_sd", "p"]].round(4).to_string(index=False))
