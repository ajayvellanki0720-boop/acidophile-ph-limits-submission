"""Why is exposed CYTOPLASMIC also positive in sister pairs? -> DATA_DIR/sister_pairs_residualised.csv
In 45 the sister-pair test on the wide panel found the exposed periplasmic class positive (q 0.0007) but also the exposed cytoplasmic class (q 0.077), while in the phylogenetic
regressions the cytoplasmic class is reliably null. The likely reason is a difference in what the two tests control: every regression here conditions on genome-wide acidic
excess, and the sister-pair test in 45 did not, so between two close relatives a shift in whole-proteome acidity would move BOTH exposed classes together.
This repeats 45 with each class's composite first residualised on genome-wide acidic excess (ordinary least squares across all bacteria, residual taken before the pair
differences), which is the same control the regressions apply. If the explanation is right, the periplasmic class should survive and the cytoplasmic class should not.
Usage: python3 scripts/46_sister_pairs_residualised.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from Bio import Phylo
from scipy import stats
warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
B = df[df.domain == "Bacteria"].drop_duplicates("ncbi_assembly").set_index("ncbi_assembly")
B["gw"] = B.acidic_excess * 100
tree = Phylo.read(T / "Bacteria.nwk", "newick")
z = lambda s: (s - s.mean()) / s.std()
pairs = []
for cl in tree.get_nonterminals():
    c = cl.clades
    if len(c) == 2 and all(x.is_terminal() for x in c):
        a, b = c[0].name, c[1].name
        if a in B.index and b in B.index and pd.notna(B.ph_min.get(a)) and pd.notna(B.ph_min.get(b)):
            hi, lo = (a, b) if B.ph_min[a] > B.ph_min[b] else (b, a)
            pairs.append((hi, lo, abs(B.ph_min[a] - B.ph_min[b])))
def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen: continue
        n = B[f"{r}|{cls}|n"]
        cols.append(z(((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
def residualise(pr):
    m = pr.notna() & B.gw.notna()
    A = np.column_stack([np.ones(int(m.sum())), B.gw[m].to_numpy(float)])
    beta, *_ = np.linalg.lstsq(A, pr[m].to_numpy(float), rcond=None)
    out = pd.Series(np.nan, index=pr.index)
    out[m] = pr[m].to_numpy(float) - A @ beta
    return out
rows = []
print(f"{len(pairs)} sister pairs; each composite residualised on genome-wide acidic excess before differencing\n")
for cls, names, lab in (("peri_exp", MEM, "exposed periplasmic"), ("peri_bur", MEM, "buried periplasmic"), ("cyto_exp", MEM, "exposed cytoplasmic"),
                        ("cyto_bur", MEM, "buried cytoplasmic"), ("TM", MEM, "transmembrane"), ("surface", SOL, "soluble surface")):
    raw = composite(names, cls)
    for tag, pr in (("raw (as in 45)", raw), ("composition removed", residualise(raw))):
        sd = pr.std()
        for md in (0.0, 0.5, 1.0):
            d = [(pr.get(hi), pr.get(lo)) for hi, lo, dd in pairs if dd > md]
            d = [(x, y) for x, y in d if pd.notna(x) and pd.notna(y)]
            if len(d) < 30: continue
            diff = np.array([x - y for x, y in d]) / sd
            pos = int((diff > 0).sum()); nz = int((diff != 0).sum())
            w = stats.wilcoxon(diff, alternative="greater") if nz > 10 else None
            rows.append({"class": lab, "control": tag, "min_dpH": md, "pairs": len(d), "mean_diff_sd": diff.mean(),
                         "frac_positive": pos / nz, "p_signrank": w.pvalue if w else np.nan})
R = pd.DataFrame(rows)
for tag, g in R.groupby("control"):
    R.loc[g.index, "q"] = stats.false_discovery_control(g.p_signrank.fillna(1).to_numpy(), method="bh")
R.to_csv(D / "sister_pairs_residualised.csv", index=False)
pd.set_option("display.width", 200)
for tag in ("raw (as in 45)", "composition removed"):
    print(f"=== {tag}")
    print(R[R.control == tag][["class", "min_dpH", "pairs", "mean_diff_sd", "frac_positive", "p_signrank", "q"]].round(4).to_string(index=False))
    print(f"    surviving q<0.1: {', '.join(sorted(set(R[(R.control == tag) & (R.q < 0.1)]['class']))) or 'none'}\n")
