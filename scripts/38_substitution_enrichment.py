"""Substitution-direction test in sister pairs -> DATA_DIR/substitution_enrichment.csv, DATA_DIR/substitution_residue_level.csv and DATA_DIR/substitution_calibration.csv.
Replaces the lab-evolution enrichment test specified in docs/CHARGE_REVERSAL_DESIGN.md section 10, which is not runnable: ALEdb exposes no API and documents no acid or low-pH
experiment, and the published acid-adaptation experiments in Escherichia coli report only a few missense mutations per clone, concentrated in regulators (rpoC, rpoB, evgS, cadA)
rather than in the conserved membrane proteins tested here, so the compartment contrast has no power there. This test asks the same question of natural evolution's own recent
substitutions instead.

Design: take every sister pair (cherry) in the bacterial marker tree - 844 pairs, 600 differing in pH_min, most within a genus - so each comparison is between closest
relatives and the difference between them reflects recent change rather than standing correlation across the tree.
  Part 1 (compartment level): within each pair, is the acidic excess of a structural class higher in the member with the HIGHER pH_min (the less acid-tolerant one)? Signed-rank
  and sign tests on the paired differences, for the exposed periplasmic class and for the buried periplasmic, exposed cytoplasmic, transmembrane and soluble-surface controls.
  Part 2 (residue level, the direct analogue of the lab-evolution test): among individual positions where the two sisters carry different residues and exactly one is acidic,
  how often does the acidic residue sit in the less acid-tolerant member? Binomial test against 50:50, per structural class. This counts natural substitutions and asks whether
  they are polarised in the predicted direction and whether that polarisation is specific to exposed periplasmic positions.
Usage: python3 scripts/38_substitution_enrichment.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
SOL = ["EF-Tu", "GroEL", "RecA"]
CLASSES = ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"]
MINDIFF = [0.0, 0.5, 1.0]
reflen = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size()
CLS = pd.read_csv(ROOT / "results/structure_panel_classes.csv")
ph = pd.read_csv(D / "phenotypes_all.csv")
pan = pd.read_csv(D / "structure_panel.csv")
df = ph.merge(pan, on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick")
B = df[df.domain == "Bacteria"].drop_duplicates("ncbi_assembly").set_index("ncbi_assembly")
pairs = []
for cl in tree.get_nonterminals():
    c = cl.clades
    if len(c) == 2 and all(x.is_terminal() for x in c):
        a, b = c[0].name, c[1].name
        if a in B.index and b in B.index and pd.notna(B.ph_min.get(a)) and pd.notna(B.ph_min.get(b)):
            hi, lo = (a, b) if B.ph_min[a] > B.ph_min[b] else (b, a)   # hi = LESS acid tolerant
            pairs.append((hi, lo, abs(B.ph_min[a] - B.ph_min[b])))
print(f"{len(pairs)} sister pairs; {sum(1 for p in pairs if p[2] > 0)} differ in pH_min", flush=True)


def excess(ref, cls):
    n = B[f"{ref}|{cls}|n"]
    need = max(8, 0.5 * reflen[(ref, cls)]) if (ref, cls) in reflen else 8
    return ((B[f"{ref}|{cls}|D"] + B[f"{ref}|{cls}|E"] - B[f"{ref}|{cls}|K"] - B[f"{ref}|{cls}|R"]) / n).where(n >= need) * 100


# ---------- Part 1: compartment level ----------
rows = []
for cls in CLASSES + ["surface"]:
    refs = SOL if cls == "surface" else MEM
    cols = [excess(r, cls) for r in refs if f"{r}|{cls}|n" in B]
    if not cols:
        continue
    M = pd.concat([(c - c.mean()) / c.std() for c in cols], axis=1)
    comp = M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
    for md in MINDIFF:
        d = [(comp.get(hi), comp.get(lo)) for hi, lo, dd in pairs if dd > md]
        d = [(x, y) for x, y in d if pd.notna(x) and pd.notna(y)]
        if len(d) < 30:
            continue
        diff = np.array([x - y for x, y in d])
        pos = int((diff > 0).sum()); nz = int((diff != 0).sum())
        w = stats.wilcoxon(diff, alternative="greater") if nz > 10 else None
        rows.append({"level": "compartment", "class": "soluble surface" if cls == "surface" else cls, "min_dpH": md, "pairs": len(d),
                     "mean_diff_sd": diff.mean(), "median_diff_sd": float(np.median(diff)), "frac_positive": pos / nz if nz else np.nan,
                     "p_signrank": w.pvalue if w else np.nan, "p_sign": stats.binomtest(pos, nz, 0.5, alternative="greater").pvalue if nz else np.nan})
R1 = pd.DataFrame(rows)
R1["q_signrank"] = stats.false_discovery_control(R1.p_signrank.fillna(1).to_numpy(), method="bh")   # 18 compartment tests
pd.set_option("display.width", 220)
print("\nPart 1 - compartment level: is the class's acidic excess higher in the LESS acid-tolerant sister?")
print("(mean_diff_sd in SD units of that class; frac_positive = share of pairs in the predicted direction; one-sided tests)")
print(R1.round(4).to_string(index=False))
print(f"  after Benjamini-Hochberg over the {len(R1)} compartment tests: {int((R1.q_signrank < 0.1).sum())} with q < 0.1")

# ---------- Part 2: residue level ----------
rows2 = []
for ref in ["YidC", "NuoL", "NuoN"]:
    f = D / f"position_residues_{ref}.csv.gz"
    if not f.exists():
        continue
    Mx = pd.read_csv(f, index_col=0)
    Mx.columns = range(Mx.shape[1])
    cl = CLS[CLS.ref == ref].sort_values("pos")
    for cls in CLASSES:
        idx = (cl[cl["class"] == cls].pos.to_numpy() - 1)
        if len(idx) == 0:
            continue
        for md in MINDIFF:
            n_hi = n_lo = 0
            for hi, lo, dd in pairs:
                if dd <= md or hi not in Mx.index or lo not in Mx.index:
                    continue
                a = Mx.loc[hi].to_numpy()[idx]; b = Mx.loc[lo].to_numpy()[idx]
                ok = (a != "-") & (b != "-") & (a != b)
                aa, ba = np.isin(a, ["D", "E"]), np.isin(b, ["D", "E"])
                n_hi += int((ok & aa & ~ba).sum())   # acidic in the less acid-tolerant member
                n_lo += int((ok & ba & ~aa).sum())   # acidic in the more acid-tolerant member
            tot = n_hi + n_lo
            if tot < 40:
                continue
            rows2.append({"level": "residue", "ref": ref, "class": cls, "min_dpH": md, "substitutions": tot,
                          "acidic_in_less_tolerant": n_hi, "acidic_in_more_tolerant": n_lo, "frac": n_hi / tot,
                          "p": stats.binomtest(n_hi, tot, 0.5, alternative="greater").pvalue})
R2 = pd.DataFrame(rows2)
print("\nPart 2 - residue level: among positions where the sisters differ and exactly one residue is acidic,")
print("how often is the acidic one in the LESS acid-tolerant member? (0.5 = no polarisation; one-sided binomial)")
print(R2.round(4).to_string(index=False))
if len(R2):
    print("\npooled over the three references:")
    for cls in CLASSES:
        z = R2[(R2["class"] == cls) & (R2.min_dpH == 0.5)]
        if len(z):
            h, l = int(z.acidic_in_less_tolerant.sum()), int(z.acidic_in_more_tolerant.sum())
            print(f"  {cls:9s} {h+l:6d} substitutions, {h/(h+l):.3f} in the predicted direction, p {stats.binomtest(h, h+l, 0.5, alternative='greater').pvalue:.2e}")
# ---------- Part 3: power calibration ----------
# A slope of pH on composition does NOT imply a large composition difference for a given pH difference. The expected difference in
# sister pairs is the REVERSE regression coefficient, slope * var(dx) / var(dy), which is small when the predictor explains little variance.
cols = [excess(r, "peri_exp") for r in MEM if f"{r}|peri_exp|n" in B]
Mp = pd.concat([(c - c.mean()) / c.std() for c in cols], axis=1)
cmp_ = Mp.mean(axis=1).where(Mp.notna().sum(axis=1) >= 3)
dd = np.array([[B.ph_min[hi] - B.ph_min[lo], cmp_[hi] - cmp_[lo]] for hi, lo, _ in pairs if pd.notna(cmp_.get(hi)) and pd.notna(cmp_.get(lo))])
dy, dx = dd[:, 0], dd[:, 1]
print("\nPart 3 - power calibration for the exposed periplasmic composite")
# pairs are oriented by pH, so dy is non-negative and dy.var() is NOT the variance of the signed difference;
# the signed differences are symmetric with mean 0, so use the orientation-invariant mean of squares.
vy, vx = float((dy ** 2).mean()), float((dx ** 2).mean())
print(f"  {len(dd)} pairs | E[dpH^2] {vy:.3f} | E[dcomposite^2] {vx:.4f}")
rows3 = []
for slope, nm in ((0.101, "between-species"), (0.312, "within-family")):
    brev = slope * vx / vy
    print(f"  {nm} slope {slope} pH/SD -> expected {brev:.4f} SD of composite per pH unit of difference")
    for md in MINDIFF:
        m = np.abs(dy) > md
        obs = (np.sign(dy[m]) * dx[m]).mean(); se = (np.sign(dy[m]) * dx[m]).std() / np.sqrt(m.sum()); exp = brev * np.abs(dy[m]).mean()
        print(f"     |dpH|>{md}: n={m.sum():3d} observed {obs:+.4f} +/- {se:.4f} SD, expected {exp:+.4f} -> {abs(obs - exp) / se:.1f} se from expectation")
        rows3.append({"level": "calibration", "slope_source": nm, "min_dpH": md, "pairs": int(m.sum()), "observed_sd": obs, "se": se, "expected_sd": exp})
pd.DataFrame(rows3).to_csv(D / "substitution_calibration.csv", index=False)
R1.to_csv(D / "substitution_enrichment.csv", index=False)
R2.to_csv(D / "substitution_residue_level.csv", index=False)
