"""Diagnostics for the broadened-sample results: are acidic excess / mean pI (new top features) and histidine (collapsed) acid-tolerance signals or an artefact of the
wider pH range (alkaliphiles, halophiles) or of GC / genome size? Subsets are chosen on PREDICTORS (reported optimum, phylum, source), not on the outcome.
Usage: python3 scripts/30_broadened_diagnostics.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["log_genome"] = np.log10(B.genome_size_bp)
C_ALL = ph.vcv_fast(tree, list(B.ncbi_assembly)); Y = B.ph_min.to_numpy(float)
FE = ["acidic_excess", "mean_pI", "frac_His", "hm_rate_kdp", "frac_acidic", "frac_basic"]
rows = []


def run(label, mask, extra=(), with_opt=False):
    idx = np.where(np.asarray(mask))[0]
    g = B.iloc[idx]
    cols = (["ph_opt"] if with_opt else []) + [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1] + list(extra)
    X = np.column_stack([np.ones(len(idx))] + [g[c].to_numpy(float) for c in cols])
    b, se, p, lam = ph.pgls_scan(C_ALL[np.ix_(idx, idx)], Y[idx], X, g[FE].to_numpy(float))
    for f, bb, ss, pp in zip(FE, b, se, p):
        rows.append({"analysis": label, "n": len(idx), "feature": f, "beta_per_sd": bb, "se": ss, "p": pp, "lambda": lam})


run("all bacteria (total effect)", np.ones(len(B), bool))
run("+ GC + genome size", np.ones(len(B), bool), ("gc", "log_genome"))
has = B.ph_opt.notna().to_numpy()
run("reported optimum <= 7.5 (alkaliphiles excluded), total effect", has & (B.ph_opt <= 7.5).to_numpy())
run("reported optimum <= 7.5, + GC + genome size", has & (B.ph_opt <= 7.5).to_numpy(), ("gc", "log_genome"))
run("reported optimum <= 7.5, optimum as covariate", has & (B.ph_opt <= 7.5).to_numpy(), with_opt=True)
run("reported optimum <= 7.5, optimum + GC + genome size", has & (B.ph_opt <= 7.5).to_numpy(), ("gc", "log_genome"), with_opt=True)
run("reported optimum <= 6.5 (acidophiles and acid-tolerant)", has & (B.ph_opt <= 6.5).to_numpy())
new = (B.origin == "literature_auto").to_numpy()
run("the earlier 1,500-species snapshot only", np.isin(B.ncbi_assembly, pd.read_csv(ROOT / "results/broad/phenotypes_all.csv").ncbi_assembly))
run("species added in this expansion only", ~np.isin(B.ncbi_assembly, pd.read_csv(ROOT / "results/broad/phenotypes_all.csv").ncbi_assembly))
R = pd.DataFrame(rows); R.to_csv(D / "broadened_diagnostics.csv", index=False)
pd.set_option("display.width", 200)
t = R.pivot_table(index="analysis", columns="feature", values="beta_per_sd", sort=False).round(3)
n = R.drop_duplicates("analysis").set_index("analysis").n
p = R.pivot_table(index="analysis", columns="feature", values="p", sort=False)
print("beta per SD (n):"); print(t.assign(n=n).to_string())
print("\np-values:"); print(p.applymap(lambda x: f"{x:.1e}").to_string())
print("\nspearman(pH_min, feature) by predictor subset for context; corr(GC, mean_pI)=%.2f corr(GC, acidic_excess)=%.2f corr(pH_min, GC)=%.2f" % (B.gc.corr(B.mean_pI), B.gc.corr(B.acidic_excess), B.ph_min.corr(B.gc)))
print("pH_min by optimum class (species with an optimum):", B[has].groupby(pd.cut(B[has].ph_opt, [0, 5, 6.5, 7.5, 8.5, 12])).ph_min.agg(["count", "mean"]).round(2).to_dict("index"))
