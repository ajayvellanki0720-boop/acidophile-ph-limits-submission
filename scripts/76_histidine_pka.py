"""Does the residue ranking on the exposed periplasmic surface follow side-chain pKa? -> DATA_DIR/histidine_pka.csv, DATA_DIR/panel_reproducibility.csv
WRITTEN BEFORE THE DATA EXIST. The panel rebuild that adds histidine takes hours, so this script and its predictions were written while it ran. That ordering is
deliberate: four times in four days a docstring here asserted a RESULT before the numbers were in and had to be weakened. Stating a PREDICTION beforehand is the opposite
move, and it is the only way this test means anything, because both outcomes are publishable and it would be trivially easy to narrate whichever one arrived.
THE QUESTION. scripts/75 found glutamate carries about 2.4x what aspartate does on the exposed periplasmic surface, and that the difference survives controlling GC. Glutamate's
side-chain pKa is about 4.25 and aspartate's about 3.65, so the stronger residue is the one whose charge state turns over at the HIGHER pH, i.e. earlier as the envelope
acidifies. That is what a charge-burden account predicts. But it is one comparison with two levels - about one degree of freedom - and pKa is confounded there with side-chain
length and flexibility. Histidine (pKa ~6.0) is the discriminating case.
THE TWO PREDICTIONS, FIXED HERE.
  A. CHARGE BURDEN. If what matters is how much a surface's charge state changes over the growth range, histidine turns over earliest of the three and should behave like a
     STRONGER GLUTAMATE: positive coefficient (more histidine -> higher pH_min, i.e. less acid tolerance), and larger than glutamate's. Ordering H > E > D.
  B. BUFFERING. Histidine is the main protein buffer precisely because its pKa sits in the physiological range, so more of it on the envelope-facing surface could BUFFER the
     periplasm and help. That predicts a NEGATIVE coefficient - the opposite sign from both carboxylates - and it is the direction this project's own withdrawn histidine lead
     ran proteome-wide (-0.55 per point across four phyla, docs/SYNTHESIS.md section 4).
  Either result discriminates, which is rare here. A null is the uninformative outcome and would most likely mean histidine is too rare on these surfaces to estimate.
WHAT WOULD MAKE THIS UNINTERPRETABLE, stated in advance. Histidine is scarce: if the exposed periplasmic class carries so few histidines that the share is mostly zero, the
coefficient is estimated off a handful of genomes with any at all and should not be read. Part 2 reports the share's distribution and refuses the test if the median reference
has under 2 histidines in that class, rather than reporting a number that looks like an answer.
PART 1 IS A REPRODUCIBILITY AUDIT AND IS WORTH MORE THAN THE REST. The rebuild writes a NEW file rather than overwriting the panel every result in this project rests on, so
the D, E, K, R and n columns of the two files can be compared directly. They should be identical: same references, same genomes, same deterministic alignment. If they are,
that is the first end-to-end reproduction of the panel that carries the headline. If they are not, every figure downstream moves, and that matters far more than histidine.
Usage: python3 scripts/76_histidine_pka.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
NEW, OLD = D / "structure_panel_wide_h.csv", D / "structure_panel_wide.csv"
MIN_MEDIAN_H = 2          # refuse the test below this, declared before seeing the data
rows = []
for f in (NEW, OLD):
    if not f.exists():
        print(f"{f.name} missing - the panel rebuild has not finished; nothing to do yet")
        sys.exit(0)
W, O = pd.read_csv(NEW), pd.read_csv(OLD)

print("1. reproducibility of the panel that carries the headline (the more important half of this script)")
shared_rows = sorted(set(W.ncbi_assembly) & set(O.ncbi_assembly))
a = W.set_index("ncbi_assembly").loc[shared_rows]
b = O.set_index("ncbi_assembly").loc[shared_rows]
cols = [c for c in b.columns if "|" in c and c.split("|")[2] in ("n", "D", "E", "K", "R")]
cols = [c for c in cols if c in a.columns]
print(f"   {len(shared_rows)} genomes in both files ({len(W)} new, {len(O)} old); {len(cols)} D/E/K/R/n columns in common")
exact = same = 0
worst = []
for c in cols:
    x, y = a[c], b[c]
    both = x.notna() & y.notna()
    eq = int((x[both] == y[both]).sum())
    same += eq
    exact += int(both.sum())
    if both.sum() and eq < both.sum():
        worst.append((c, int(both.sum()) - eq, int(both.sum())))
frac = same / exact if exact else float("nan")
print(f"   {same} of {exact} cell comparisons identical ({frac:.6f})")
if worst:
    print(f"   {len(worst)} columns differ somewhere; the five largest discrepancies:")
    for c, nd, nt in sorted(worst, key=lambda z: -z[1])[:5]:
        print(f"      {c:34s} {nd} of {nt} cells differ")
    print("   THIS IS THE HEADLINE OF THIS SCRIPT, not histidine: the panel does not rebuild identically, so the")
    print("   figures every result rests on depend on something not captured by the inputs. Investigate before")
    print("   reading part 3 at all.")
else:
    print("   Every shared cell is identical. This is the first end-to-end reproduction of the 32-reference panel,")
    print("   and it means the rebuild adds histidine without disturbing anything the project already rests on.")
missing_cols = [c for c in O.columns if "|" in c and c.split("|")[2] in ("n", "D", "E", "K", "R") and c not in W.columns]
if missing_cols:
    print(f"   NOTE: {len(missing_cols)} old columns are absent from the new file (e.g. {missing_cols[:3]})")
rows.append({"part": "reproducibility", "genomes_shared": len(shared_rows), "genomes_new": len(W),
             "genomes_old": len(O), "columns_compared": len(cols), "cells_compared": exact,
             "cells_identical": same, "fraction_identical": frac, "columns_differing": len(worst),
             "old_columns_absent": len(missing_cols)})
pd.DataFrame(rows).to_csv(D / "panel_reproducibility.csv", index=False)

print("\n2. is histidine common enough on these surfaces to estimate at all?")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in refs[refs.accepted & refs.membrane].label]
hcols = [f"{r}|peri_exp|H" for r in MEM if f"{r}|peri_exp|H" in W.columns]
ncols = [f"{r}|peri_exp|n" for r in MEM if f"{r}|peri_exp|H" in W.columns]
med_h = float(np.nanmedian([W[c].median() for c in hcols])) if hcols else 0.0
print(f"   {len(hcols)} references carry a histidine column for the exposed periplasmic class")
if hcols:
    print(f"   median histidines per genome, median over references: {med_h:.2f}")
    print(f"   median class size (aligned positions): {np.nanmedian([W[c].median() for c in ncols]):.1f}")
rows.append({"part": "histidine abundance", "references": len(hcols), "median_H_per_genome": med_h,
             "threshold": MIN_MEDIAN_H})
if not hcols or med_h < MIN_MEDIAN_H:
    print(f"   BELOW THE PRE-DECLARED FLOOR of {MIN_MEDIAN_H}: histidine is too scarce on these surfaces for a share to")
    print("   be estimated, so the test is REFUSED rather than reported. This was declared before the data existed;")
    print("   see the docstring. The pKa question is not answered either way and needs a different surface or panel.")
    pd.DataFrame(rows).to_csv(D / "histidine_pka.csv", index=False)
    print("\nwrote", D / "histidine_pka.csv")
    sys.exit(0)

print("\n3. histidine against the two predictions, with the carboxylates for comparison")
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp"]],
             on="ncbi_assembly")
      .merge(W, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
B["gcz"] = (B.gc - B.gc.mean()) / B.gc.std()
B["szz"] = (B.genome_size_bp - B.genome_size_bp.mean()) / B.genome_size_bp.std()
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
CLS = pd.read_csv(ROOT / "results/structure_panel_wide_h_classes.csv")
reflen = CLS.groupby(["ref", "class"]).size()


def pool(bv, sv):
    bv, sv = np.asarray(bv, float), np.asarray(sv, float)
    w = 1 / sv ** 2
    fe = (w * bv).sum() / w.sum()
    Q = (w * (bv - fe) ** 2).sum()
    k = len(bv) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (sv ** 2 + tau2)
    mu = (wr * bv).sum() / wr.sum()
    sem = float(np.sqrt(1 / wr.sum()))
    return float(mu), sem, float(2 * stats.norm.sf(abs(mu / sem)))


for cls in ("peri_exp", "cyto_exp"):
    got = []
    for ref in MEM:
        if f"{ref}|{cls}|H" not in B:
            continue
        n = B[f"{ref}|{cls}|n"]
        need = max(8, 0.5 * reflen[(ref, cls)]) if (ref, cls) in reflen.index else 8
        xs = {k: (B[f"{ref}|{cls}|{k}"] / n).where(n >= need) * 100 for k in ("D", "E", "H")}
        m = (xs["D"].notna() & xs["E"].notna() & xs["H"].notna() & B.gc.notna() & B.genome_size_bp.notna()).to_numpy()
        if m.sum() < 200 or xs["H"][m].std() == 0:
            continue
        g = B[m].reset_index(drop=True)
        Cm = Call[np.ix_(m, m)]
        keep = [k for k in SRCC if g[k].nunique() > 1]
        Xb = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw", "gcz", "szz"] + keep])
        lam = ph.pgls(g.ph_min.to_numpy(float), Xb, Cm)["lambda"]
        Xf = np.column_stack([Xb] + [xs[k][m].to_numpy() for k in ("D", "E", "H")])
        r = ph.pgls(g.ph_min.to_numpy(float), Xf, Cm, lam=lam)
        got.append({"ref": ref, "n": int(m.sum()),
                    "beta_D": float(r["beta"][-3]), "se_D": float(r["se"][-3]),
                    "beta_E": float(r["beta"][-2]), "se_E": float(r["se"][-2]),
                    "beta_H": float(r["beta"][-1]), "se_H": float(r["se"][-1])})
    if len(got) < 2:
        print(f"   {cls}: fewer than two references estimable, skipped")
        continue
    Gg = pd.DataFrame(got)
    print(f"   {cls} ({len(Gg)} references), all three residues entered together, per percentage point:")
    res = {}
    for res_name, pka in (("D", 3.65), ("E", 4.25), ("H", 6.00)):
        mu, sem, pv = pool(Gg[f"beta_{res_name}"], Gg[f"se_{res_name}"])
        res[res_name] = mu
        print(f"      {res_name} (pKa {pka:.2f})  {mu:+.5f} +/- {sem:.5f} (p {pv:.2e}); positive in "
              f"{int((Gg[f'beta_{res_name}'] > 0).sum())}/{len(Gg)}")
        rows.append({"part": "joint DEH", "class": cls, "residue": res_name, "pka": pka,
                     "beta_per_point": mu, "se": sem, "p": pv, "n_refs": int(len(Gg)),
                     "n_refs_positive": int((Gg[f"beta_{res_name}"] > 0).sum())})
    if cls == "peri_exp":
        ordered = res["H"] > res["E"] > res["D"] > 0
        opposite = res["H"] < 0 < res["E"]
        verdict = ("A: charge burden - ordering follows pKa" if ordered else
                   "B: buffering - histidine takes the opposite sign" if opposite else
                   "NEITHER prediction cleanly: see the coefficients")
        print(f"      VERDICT, against the predictions fixed in the docstring: {verdict}")
        rows.append({"part": "verdict", "class": cls, "beta_D": res["D"], "beta_E": res["E"], "beta_H": res["H"],
                     "prediction_A_charge_burden": bool(ordered), "prediction_B_buffering": bool(opposite),
                     "verdict": verdict})
print("   The cytoplasmic class is the specificity control: a pKa effect should be weaker or absent on the side")
print("   that does not face the acidified medium, and anything that moves both classes alike is not about pH.")

pd.DataFrame(rows).to_csv(D / "histidine_pka.csv", index=False)
print("\nwrote", D / "histidine_pka.csv", "and", D / "panel_reproducibility.csv")
print("\nReading. Part 1 governs: if the panel does not rebuild cell-for-cell, nothing else here is worth reading and")
print("every downstream figure needs re-examining. If it does, part 3 is a genuine test between two mechanisms rather")
print("than a confirmation of either, because the predictions were fixed before the data existed and both were")
print("publishable. A null means histidine is too scarce on these surfaces, which part 2 refuses in advance rather")
print("than dressing up as an answer.")
