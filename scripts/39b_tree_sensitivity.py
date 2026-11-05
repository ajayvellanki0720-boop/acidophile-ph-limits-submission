"""Is the exposed-periplasmic result sensitive to the tree? -> DATA_DIR/tree_sensitivity.csv and DATA_DIR/tree_sensitivity_wide.csv
Every result on the 2,768-genome set so far rests on one neighbour-joining tree built from the concatenated ribosomal-protein alignment, with no replicates. This re-fits the key
PGLS on 25 column-bootstrap trees of the same alignment (built by scripts/08 into results/phylo_broad2_boot, leaving the committed tree untouched) and reports the spread of the
estimate across them against the point estimate on the main tree. Run for the exposed periplasmic composite and, as the matched control, the soluble-surface composite, each
with genome-wide acidic excess controlled. A result that moves little across bootstrap trees is not an artefact of one topology; it does not address systematic error shared by
all neighbour-joining trees of this alignment, which would need a maximum-likelihood tree.
Usage: python3 scripts/39b_tree_sensitivity.py [panel]   panel = 'orig' (default) or 'wide'
       (env DATA_DIR=results/broad2, BOOT_DIR=results/phylo_broad2_boot, PHYLO_DIR=results/phylo_broad2)"""
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
WHICH = sys.argv[1] if len(sys.argv) > 1 else "orig"
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
BOOT = ROOT / os.environ.get("BOOT_DIR", "results/phylo_broad2_boot")
MAIN = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
if WHICH == "wide":
    refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
    MEM = [str(r.label) for r in refs.itertuples() if r.membrane]; SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
    panel_csv, classes = D / "structure_panel_wide.csv", ROOT / "results/structure_panel_wide_classes.csv"
else:
    MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]; SOL = ["EF-Tu", "GroEL", "RecA"]
    panel_csv, classes = D / "structure_panel.csv", ROOT / "results/structure_panel_classes.csv"
reflen = pd.read_csv(classes).groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(panel_csv), on="ncbi_assembly"))
main = Phylo.read(MAIN / "Bacteria.nwk", "newick"); tips = {t.name for t in main.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]


def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen:
            continue
        n = B[f"{r}|{cls}|n"]
        x = ((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))
        cols.append((x - x.mean()) / x.std())
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)


PRED = {"exposed periplasmic": composite(MEM, "peri_exp"), "soluble surface": composite(SOL, "surface")}
boots = [t for t in Phylo.parse(str(BOOT / "Bacteria_bootstrap.nwk"), "newick")]
print(f"panel '{WHICH}': {len(MEM)} membrane, {len(SOL)} soluble | {len(boots)} bootstrap trees | {len(B)} bacteria\n", flush=True)
rows = []
for name, pr in PRED.items():
    m = pr.notna().to_numpy()
    g = B[m].reset_index(drop=True); keep = list(g.ncbi_assembly)
    # pgls_scan standardises the tested predictor internally, so its beta is ALREADY per SD; do not rescale it.
    X = np.column_stack([np.ones(len(g)), g.gw] + [g[c] for c in SRCC])
    y = g.ph_min.to_numpy(float); xv = pr[m].to_numpy(float)[:, None]
    C0 = ph.vcv_fast(main, keep)
    b0, se0, p0, lam0 = ph.pgls_scan(C0, y, X, xv)
    print(f"{name}: main tree {b0[0]:+.4f} per SD (se {se0[0]:.4f}, p {p0[0]:.2e}, lambda {lam0:.2f}), n={len(g)}", flush=True)
    est = []
    for i, t in enumerate(boots):
        try:
            Ci = ph.vcv_fast(t, keep)
        except Exception as e:
            print(f"   tree {i}: skipped ({str(e)[:40]})", flush=True); continue
        bi, sei, pi, lami = ph.pgls_scan(Ci, y, X, xv)
        est.append((bi[0], sei[0], pi[0], lami))
        rows.append({"panel": WHICH, "predictor": name, "tree": i, "beta_per_sd": bi[0], "se": sei[0], "p": pi[0], "lambda": lami})
    e = np.array(est)
    print(f"   across {len(e)} bootstrap trees: median {np.median(e[:, 0]):+.4f}, range {e[:, 0].min():+.4f} to {e[:, 0].max():+.4f}; "
          f"p < 0.05 in {int((e[:, 2] < 0.05).sum())}/{len(e)}; lambda {np.median(e[:, 3]):.2f}", flush=True)
    rows.append({"panel": WHICH, "predictor": name, "tree": -1, "beta_per_sd": b0[0], "se": se0[0], "p": p0[0], "lambda": lam0})
out = D / f"tree_sensitivity{'_wide' if WHICH == 'wide' else ''}.csv"
# the median and range across bootstrap trees are what the synthesis quotes, and they were printed but never
# stored, so they were untraceable to this file. A derived figure quoted anywhere has to be a stored row.
_R = pd.DataFrame(rows)
_summ = []
for _pr, _g in _R[_R.tree >= 0].groupby("predictor"):
    # min and max get their OWN columns. Overloading se and p to carry them would plant exactly the kind of
    # mislabelled-units trap this project has had to correct five times.
    _summ.append({"panel": WHICH, "predictor": _pr, "tree": -2, "beta_per_sd": float(_g.beta_per_sd.median()),
                  "beta_median_across_trees": float(_g.beta_per_sd.median()),
                  "beta_min_across_trees": float(_g.beta_per_sd.min()),
                  "beta_max_across_trees": float(_g.beta_per_sd.max()),
                  "n_trees": int(len(_g)), "n_significant": int((_g.p < 0.05).sum()),
                  "lambda": float(_g["lambda"].median()),
                  "note": "tree=-2 is a summary row across the bootstrap trees, not a fit"})
pd.concat([_R, pd.DataFrame(_summ)], ignore_index=True).to_csv(out, index=False)
print("\nwrote", out)
