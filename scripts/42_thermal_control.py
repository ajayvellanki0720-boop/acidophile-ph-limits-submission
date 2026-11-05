"""Is the exposed-periplasmic result a thermophily artefact? -> DATA_DIR/thermal_control.csv
Thermophily was never controlled in the exposure tests, and it is a plausible confound: thermoacidophiles are both hot and acid-tolerant, growth temperature correlates -0.16
with pH_min across these bacteria, and growth temperature reshapes proteome composition. Two thermal covariates are used because each has a weakness. Growth temperature
optimum is the real variable but is known for only 1,467 of 2,771 bacteria, mostly from TEMPURA. IVYWREL, the proteome fraction of the residues that track thermal adaptation,
is available for every genome and correlates 0.79 with the optimum here, but it is itself a composition measure, so controlling for it is conservative to the point of removing
some of the signal of interest.
Tests on the wide 32-protein panel composite and its matched soluble-surface control: thermal covariates added singly and together, then restricted to mesophiles
(optimum below 45 C) and to thermophiles separately, so the association can be seen inside a temperature band rather than across one.
Usage: python3 scripts/42_thermal_control.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
feat = pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp", "ivywrel"]]
df = pd.read_csv(D / "phenotypes_all.csv").merge(feat, on="ncbi_assembly").merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
z = lambda s: (s - s.mean()) / s.std()
B["gw"] = B.acidic_excess * 100
B["gcz"], B["sz"], B["ivy"] = z(B.gc), z(B.genome_size_bp), z(B.ivywrel)
B["topt"] = z(B.temp_opt_c)
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
print(f"{len(B)} bacteria | growth optimum temperature known for {int(B.temp_opt_c.notna().sum())} | "
      f"corr(IVYWREL, temp_opt) {B[['ivywrel', 'temp_opt_c']].corr().iloc[0, 1]:.2f}\n")


def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen:
            continue
        n = B[f"{r}|{cls}|n"]
        x = ((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))
        cols.append(z(x))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)


PRED = {"exposed periplasmic": composite(MEM, "peri_exp"), "soluble surface": composite(SOL, "surface")}
rows = []
MODELS = [("+ composition (reference)", ["gw"], None),
          ("+ composition, IVYWREL", ["gw", "ivy"], None),
          ("+ composition, growth temperature", ["gw", "topt"], "topt"),
          ("+ composition, IVYWREL, growth temperature", ["gw", "ivy", "topt"], "topt"),
          ("+ composition, GC, size, IVYWREL, temperature", ["gw", "gcz", "sz", "ivy", "topt"], "topt")]
SUBSETS = [("all bacteria", None), ("mesophiles only (optimum < 45 C)", lambda d: d.temp_opt_c < 45), ("thermophiles only (optimum >= 45 C)", lambda d: d.temp_opt_c >= 45)]
for name, pr in PRED.items():
    for sub, f in SUBSETS:
        for lab, cov, need in MODELS:
            if sub != "all bacteria" and need is None and lab != "+ composition (reference)":
                continue
            m = pr.notna().to_numpy()
            if need is not None:
                m &= B.temp_opt_c.notna().to_numpy()
            if f is not None:
                m &= f(B).fillna(False).to_numpy()
            if m.sum() < 120:
                continue
            g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
            use = [c for c in cov + SRCC if g[c].nunique() > 1]
            X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
            b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, pr[m].to_numpy(float)[:, None])
            rows.append({"predictor": name, "subset": sub, "model": lab, "n": len(g), "beta_per_sd": b[0], "se": se[0], "p": p[0]})
R = pd.DataFrame(rows); R.to_csv(D / "thermal_control.csv", index=False)
pd.set_option("display.width", 220)
for name in PRED:
    print(f"=== {name} (beta per SD of the composite)")
    print(R[R.predictor == name][["subset", "model", "n", "beta_per_sd", "se", "p"]].round(4).to_string(index=False))
    print()
