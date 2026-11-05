"""How much does the composite actually add to a prediction of pH_min? -> DATA_DIR/predictive_value.csv
Everything so far is inference about a coefficient. This asks the practical question instead: out of sample, does knowing the exposed-periplasmic composite improve a prediction
of the lower limit over knowing the growth optimum alone? The earlier full-dataset models concluded that no genomic model beats the optimum (R-squared 0.82), but that was
before this predictor existed.
Cross-validation is blocked by FAMILY - whole families are held out together - because species inside a family are not independent and a random split would let the model
memorise family-level means and report a flattering number. Ten folds, families assigned to folds at random, with the composite standardised on the training fold only. Models:
the optimum alone, the composite alone, both together, and both plus genome-wide acidic excess. Reported as out-of-sample R-squared against the test-fold mean, so a negative
value means worse than predicting the average.
Usage: python3 scripts/49_predictive_value.py [n_repeats]   (env DATA_DIR=results/broad2)"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
B = df[df.domain == "Bacteria"].drop_duplicates("ncbi_assembly").reset_index(drop=True)
z = lambda s: (s - s.mean()) / s.std()
cols = []
for r in MEM:
    if f"{r}|peri_exp|n" not in B or (r, "peri_exp") not in reflen: continue
    n = B[f"{r}|peri_exp|n"]
    cols.append(z(((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, "peri_exp")]))))
M = pd.concat(cols, axis=1)
B["x"] = M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
B["gw"] = B.acidic_excess * 100
d = B[B.x.notna() & B.ph_min.notna() & B.ph_opt.notna() & B.family.notna()].reset_index(drop=True)
print(f"{len(d)} bacteria with the composite, a minimum and an optimum, in {d.family.nunique()} families\n")
MODELS = {"optimum alone": ["ph_opt"], "composite alone": ["x"], "optimum + composite": ["ph_opt", "x"],
          "optimum + composite + composition": ["ph_opt", "x", "gw"], "optimum + composition": ["ph_opt", "gw"]}
fams = d.family.unique()
res = {k: [] for k in MODELS}
rng = np.random.default_rng(49)
for rep in range(NREP):
    assign = {f: i % 10 for i, f in enumerate(rng.permutation(fams))}
    fold = d.family.map(assign).to_numpy()
    for name, cs in MODELS.items():
        num, den = 0.0, 0.0
        for k in range(10):
            tr, te = fold != k, fold == k
            if te.sum() < 5 or tr.sum() < 50: continue
            Xtr = d.loc[tr, cs].to_numpy(float); Xte = d.loc[te, cs].to_numpy(float)
            mu, sg = Xtr.mean(0), Xtr.std(0); sg[sg == 0] = 1
            Xtr = np.column_stack([np.ones(tr.sum()), (Xtr - mu) / sg]); Xte = np.column_stack([np.ones(te.sum()), (Xte - mu) / sg])
            ytr = d.loc[tr, "ph_min"].to_numpy(float); yte = d.loc[te, "ph_min"].to_numpy(float)
            beta = np.linalg.lstsq(Xtr, ytr, rcond=None)[0]
            pred = Xte @ beta
            num += float(((yte - pred) ** 2).sum()); den += float(((yte - ytr.mean()) ** 2).sum())
        res[name].append(1 - num / den)
rows = []
print("out-of-sample R-squared, 10 folds blocked by family, averaged over", NREP, "random fold assignments")
for name in MODELS:
    v = np.array(res[name])
    rows.append({"model": name, "r2_mean": v.mean(), "r2_sd": v.std(), "n": len(d)})
    print(f"  {name:36s} {v.mean():+.4f}  (sd over repeats {v.std():.4f})")
base = np.array(res["optimum alone"]); add = np.array(res["optimum + composite"])
print(f"\nadding the composite to the optimum changes out-of-sample R-squared by {(add - base).mean():+.4f} "
      f"(sd {(add - base).std():.4f}; improves in {int((add > base).sum())}/{NREP} fold assignments)")
rows.append({"model": "delta from adding the composite", "r2_mean": float((add - base).mean()), "r2_sd": float((add - base).std()), "n": len(d)})
pd.DataFrame(rows).to_csv(D / "predictive_value.csv", index=False)
