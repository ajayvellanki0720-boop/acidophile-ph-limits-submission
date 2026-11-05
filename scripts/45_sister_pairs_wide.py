"""The sister-pair substitution-direction test, redone on the wide panel -> DATA_DIR/sister_pairs_wide.csv
Script 38 ran this on the hand-picked 7-protein composite and found nothing, concluding from its power calibration that the test was never capable of detecting the effect. Two
things have changed since: the composite now averages 32 references and so carries much less measurement error (its standardized effect is 1.8x larger), and the within-genus
test in 44 IS significant, which makes a null at the sister-pair scale worth re-examining rather than assuming.
Same design as 38: every sister pair (cherry) of the bacterial marker tree, oriented so that `hi` is the member with the HIGHER pH_min (less acid tolerant), asking whether the
class's acidic excess is higher in that member. Signed-rank and sign tests, at three thresholds on the pH difference, with the matched control classes, Benjamini-Hochberg over
all tests, and the same reverse-regression calibration so the result can be read against what the comparative slope actually predicts.
Usage: python3 scripts/45_sister_pairs_wide.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
MINDIFF = [0.0, 0.5, 1.0]
df = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly")
B = df[df.domain == "Bacteria"].drop_duplicates("ncbi_assembly").set_index("ncbi_assembly")
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
print(f"{len(pairs)} sister pairs, {sum(1 for p in pairs if p[2] > 0)} differing in pH_min\n")
def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen: continue
        n = B[f"{r}|{cls}|n"]
        cols.append(z(((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
rows = []
for cls, names, lab in (("peri_exp", MEM, "exposed periplasmic"), ("peri_bur", MEM, "buried periplasmic"), ("cyto_exp", MEM, "exposed cytoplasmic"),
                        ("cyto_bur", MEM, "buried cytoplasmic"), ("TM", MEM, "transmembrane"), ("surface", SOL, "soluble surface")):
    pr = composite(names, cls); sd = pr.std()
    for md in MINDIFF:
        d = [(pr.get(hi), pr.get(lo)) for hi, lo, dd in pairs if dd > md]
        d = [(x, y) for x, y in d if pd.notna(x) and pd.notna(y)]
        if len(d) < 30: continue
        diff = np.array([x - y for x, y in d]) / sd          # in SD units of the composite
        pos = int((diff > 0).sum()); nz = int((diff != 0).sum())
        w = stats.wilcoxon(diff, alternative="greater") if nz > 10 else None
        rows.append({"class": lab, "min_dpH": md, "pairs": len(d), "mean_diff_sd": diff.mean(), "frac_positive": pos / nz,
                     "p_signrank": w.pvalue if w else np.nan, "p_sign": stats.binomtest(pos, nz, 0.5, alternative="greater").pvalue})
R = pd.DataFrame(rows)
R["q_signrank"] = stats.false_discovery_control(R.p_signrank.fillna(1).to_numpy(), method="bh")
R.to_csv(D / "sister_pairs_wide.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))
print(f"\nsurviving Benjamini-Hochberg at q<0.1: {int((R.q_signrank < 0.1).sum())} of {len(R)} tests")
pr = composite(MEM, "peri_exp"); sd = pr.std()
dd = np.array([[B.ph_min[hi] - B.ph_min[lo], (pr.get(hi) - pr.get(lo)) / sd] for hi, lo, _ in pairs if pd.notna(pr.get(hi)) and pd.notna(pr.get(lo))])
vy, vx = float((dd[:, 0] ** 2).mean()), float((dd[:, 1] ** 2).mean())
print(f"\ncalibration (exposed periplasmic): E[dpH^2] {vy:.3f}, E[dcomposite^2] {vx:.3f} SD^2")
for slope, nm in ((0.185, "PGLS + composition, wide panel"), (0.242, "within families, raw scale")):
    brev = slope * vx / vy
    print(f"  {nm} = {slope} -> expected {brev:.4f} SD per pH unit")
    for md in MINDIFF:
        m = np.abs(dd[:, 0]) > md
        obs = (np.sign(dd[m, 0]) * dd[m, 1]).mean(); se = (np.sign(dd[m, 0]) * dd[m, 1]).std() / np.sqrt(m.sum()); exp = brev * np.abs(dd[m, 0]).mean()
        print(f"     |dpH|>{md}: n={m.sum():3d} observed {obs:+.4f} +/- {se:.4f}, expected {exp:+.4f} -> {abs(obs - exp) / se:.1f} se")
