"""Exposure test (34a/34b) on the optimum-matched subset: only genomes with a pH_opt AND every class of that protein scored, so the three models are fitted on identical genomes:
total effect, linear optimum, flexible (cubic) optimum. Also reports the fraction of the total effect remaining. Usage: python3 scripts/34d_structure_exposure_optimum_matched.py"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from Bio import Phylo
from acidophile import phylo as ph
warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
CL = {"YidC": ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"], "SecY": ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"], "EF-Tu": ["surface", "buried"]}
df = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "structure_exposure.csv"), on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_opt.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC: B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
z = (B.ph_opt - B.ph_opt.mean()) / B.ph_opt.std(); B["o1"], B["o2"], B["o3"] = z, z ** 2, z ** 3
rows = []
for ref, classes in CL.items():
    for c in classes:
        n = B[f"{ref}|{c}|n"]; B[f"{ref}_{c}"] = ((B[f"{ref}|{c}|D"] + B[f"{ref}|{c}|E"] - B[f"{ref}|{c}|K"] - B[f"{ref}|{c}|R"]) / n).where(n >= 15) * 100
    cols = [f"{ref}_{c}" for c in classes]; m = B[cols].notna().all(axis=1).to_numpy()
    g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
    src = [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
    for lab, opt in (("total effect", []), ("linear optimum", ["o1"]), ("cubic optimum", ["o1", "o2", "o3"])):
        X = np.column_stack([np.ones(len(g))] + [g[k] for k in opt + src])
        for c, col in zip(classes, cols):
            b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, g[[col]].to_numpy(float))
            rows.append({"ref": ref, "class": c, "model": lab, "n": len(g), "beta_per_pt": b[0], "se": se[0], "p": p[0]})
R = pd.DataFrame(rows); R.to_csv(D / "structure_exposure_optimum_matched.csv", index=False)
W = R.pivot_table(index=["ref", "class", "n"], columns="model", values="beta_per_pt").reset_index()
W["share_left_linear"] = W["linear optimum"] / W["total effect"]; W["share_left_cubic"] = W["cubic optimum"] / W["total effect"]
P = R.pivot_table(index=["ref", "class"], columns="model", values="p").add_prefix("p ").reset_index()
pd.set_option("display.width", 220)
print(W.merge(P, on=["ref", "class"])[["ref", "class", "n", "total effect", "linear optimum", "cubic optimum", "share_left_linear", "share_left_cubic", "p total effect", "p linear optimum", "p cubic optimum"]].round(4).to_string(index=False))
