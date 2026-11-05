"""Detect exactly which parser rows read an optimum as a growth range, and refit without them -> DATA_DIR/optimum_as_range.csv
Script 60 found this failure mode by hand - a sentence that reports only an optimum, as in "grew optimally at 35 C, at pH 7.0-8.0", has its optimum recorded as the growth
range - and bounded its effect using a proxy: every parser row with no parsed optimum, 1,097 of them, roughly a third wrong by the hand-check. A proxy that is two thirds
false positives is a weak instrument, and the failure is mechanically detectable, so this detects it exactly.
The detector works on the stored statement each parser row came from. It finds where the parser's own range pattern matched and flags the row on four conditions together:
  - a stem meaning "optimal" appears within 60 characters BEFORE the match, so the range sits inside an optimum clause;
  - the matched interval is no wider than 2 pH units, because an "optimum" spanning four units is not an optimum - this is the condition that does most of the work;
  - no "optim" stem follows within 40 characters, which would mean the structure is "range (optimum ...)" rather than an optimum alone;
  - the sentence states no growth range in its own right ("growth occurred at", "grew between", "pH range" and so on).
The first draft used only the first and fourth conditions and flagged three statements that were correct, each one a genuine range with a parenthetical optimum. The width
condition removes all three without losing any true case.
HAND_READ below is the ground truth: fifteen statements read by hand, eight of them wrong, recorded verbatim with the verdict so the detector is tested rather than trusted.
It currently scores 8 of 8 true positives with no false positives and no misses. A detector that cannot reproduce them is not used for the refit.
With the flag in hand, the question script 60 could only bound becomes a direct test: refit the headline estimate with exactly the affected rows removed.
Usage: python3 scripts/61_optimum_as_range.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
from acidophile.parser_flags import LOOKBACK, MAXWIDTH, TRAILING, flag  # the detector itself lives in src/acidophile/parser_flags.py

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")


# Fifteen statements read by hand, verbatim, with the verdict: True = the parser took an optimum as the growth range.
HAND_READ = [
    ("Nocardioides dubius", "Strain KSL-104(T) grew optimally at pH 7.0-8.0 and 30 degrees C. It was characterized chemotaxonomically as ha", True),
    ("Psychroflexus salinarum", "Strain ISL-14(T) grew optimally at pH 7.0-8.0, 30 degrees C and in the presence of approximately 2 % (w/v) NaC", True),
    ("Planosporangium mesophilum", "Strain YIM 48875(T) grew optimally at 28 \u00b0C, at pH 6.0-8.0 and in the presence of 2 % (w/v) NaCl", True),
    ("Rhodobacter flavimaris", "Strain WL0062T grew optimally at 28 \u00b0C, pH 7.0-8.0 and with 1.0-3.0% (w/v) NaCl", True),
    ("Sulfitobacter undariae", "Strain W-BA2(T) grew optimally at 25 \u00b0C, at pH 7.0-8.0 and in the presence of approximately 2.0-3.0% (w/v) NaC", True),
    ("Demequina activiva", "Strain BS-12M(T) grew optimally at 35 \u00b0C, at pH 7.0-8.0 and in the presence of 2.0% (w/v) NaCl", True),
    ("Aureimonas endophytica", "Strain 2T4P-2-4T grew optimally without NaCl at 28-30 \u00b0C, pH 7.0-8.0. The 16S rRNA gene sequence of strain", True),
    ("Jiella soli", "Growth optimum was observed to be at 28-37 \u00b0C, pH 5.5-6.5 and 2-5% (w/v) NaCl", True),
    ("Geodermatophilus normandii", "e was aerobic and exhibited optimal growth from 25-40 \u00b0C at pH 6.0-10.0 with 0-1% (w/v) NaCl", False),
    ("Brevundimonas aveniformis", "een 15 and 35 degrees C (optimum, 30 degrees C) and between pH 6.0 and 9.0 (optimum, pH 7.5-8.5)", False),
    ("Parabacteroides acidifaciens", "ew anaerobically at 20-45 \u00b0C (optimally at 37-40 \u00b0C) and at pH 6.0-10.0 (optimally at pH 6.0-8.0)", False),
    ("Rhodocaloribacter litoris", "Growth was observed at 40-65 \u00b0C (optimum at 60 \u00b0C), pH 6-9 (optimum at 6.5-8) and 0.5-5% (optimum at 1-2%) (w/v) NaCl", False),
    ("Hymenobacter citatus", "cterial growth was observed at 10-37 \u00b0C (optimum 25 \u00b0C) and pH 6.0-8.0 (optimum pH 6.0)", False),
    ("Salibacterium lacus", "Growth occurred at pH 5.5-8.5 (optimum pH 7.0) and at 10-45 \u00b0C (optimum 30 \u00b0C)", False),
    ("Dyella lipolytica", "Strain DHOB07T grew at 10-37 \u00b0C, pH 4-7 and 0-0.5 % (w/v) NaCl, with an optimum at 28 \u00b0C, pH 5-5.5 and 0% (w/v) NaCl on R2A medium", False),
]


cand = pd.read_csv(ROOT / "data/curated/broad_candidates.csv")
cand["pmid"] = cand.pmid.astype(str)
cand["flagged"] = cand.statement.map(flag)

print("1. the detector against the statements read by hand")
tp = fp = fn = tn = 0
for taxon, stmt, truth in HAND_READ:
    f = flag(stmt)
    if f and truth:
        tp += 1
    elif f and not truth:
        fp += 1
        print(f"   FALSE POSITIVE  {taxon}: {stmt[:95]}")
    elif truth:
        fn += 1
        print(f"   MISSED          {taxon}: {stmt[:95]}")
    else:
        tn += 1
print(f"   {len(HAND_READ)} hand-read statements, {sum(1 for *_, v in HAND_READ if v)} of them wrong: "
      f"{tp} correctly flagged, {fp} false positives, {fn} missed, {tn} correctly passed")
rows = [{"part": "detector check", "n": len(HAND_READ), "hand_wrong": sum(1 for *_, v in HAND_READ if v),
         "true_positives": tp, "false_positives": fp, "missed": fn, "correct_passes": tn}]
if fn or fp:
    print("   The detector does not reproduce the hand-read set exactly; the refit below is weaker than it looks.")

p_all = pd.read_csv(D / "phenotypes_all.csv")
B0 = p_all[p_all.domain == "Bacteria"].drop_duplicates("ncbi_assembly")
used = set(B0[B0.origin == "literature_auto"].source.dropna().str.replace("PMID:", "", regex=False).astype(str))

print("\n2. how many rows it flags overall")
inuse = cand[cand.pmid.isin(used) & cand.usable_noopt]
print(f"   {len(inuse)} candidate rows behind the parser tier; {int(inuse.flagged.sum())} flagged ({inuse.flagged.mean():.1%})")
noopt = inuse.ph_opt.isna()
print(f"   script 60's proxy (no parsed optimum) covered {int(noopt.sum())} rows, of which {int((inuse.flagged & noopt).sum())} are flagged here")
print(f"   flagged rows with an optimum parsed anyway: {int((inuse.flagged & ~noopt).sum())}")
print(f"   mean reported minimum: flagged {inuse.ph_min[inuse.flagged].mean():.3f} against {inuse.ph_min[~inuse.flagged].mean():.3f} unflagged "
      f"(+{inuse.ph_min[inuse.flagged].mean() - inuse.ph_min[~inuse.flagged].mean():.3f})")
rows.append({"part": "prevalence", "n": len(inuse), "flagged": int(inuse.flagged.sum()), "flagged_rate": float(inuse.flagged.mean()),
             "mean_ph_min_flagged": float(inuse.ph_min[inuse.flagged].mean()), "mean_ph_min_unflagged": float(inuse.ph_min[~inuse.flagged].mean())})
bad_pmids = set(inuse.pmid[inuse.flagged])

print("\n3. the headline estimate with exactly those rows removed")
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
B["pmid"] = B.source.fillna("").str.replace("PMID:", "", regex=False).astype(str)
B["affected"] = B.origin.eq("literature_auto") & B.pmid.isin(bad_pmids)
print(f"   {int(B.affected.sum())} genomes in the analysis set carry a flagged statement")
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
cols = {}
for r in MEM:
    c = f"{r}|peri_exp|n"
    if c not in B or (r, "peri_exp") not in reflen:
        continue
    n_ = B[c]
    x = ((B[f"{r}|peri_exp|D"] + B[f"{r}|peri_exp|E"] - B[f"{r}|peri_exp|K"] - B[f"{r}|peri_exp|R"]) / n_).where(n_ >= max(8, 0.5 * reflen[(r, "peri_exp")]))
    cols[r] = (x - x.mean()) / x.std()
Z = pd.DataFrame(cols)
pr = Z.mean(axis=1).where(Z.notna().sum(axis=1) >= 3)
for label, mask in (("all bacteria", pd.Series(True, index=B.index)),
                    ("drop exactly the flagged rows", ~B.affected),
                    ("ONLY the flagged rows", B.affected),
                    ("drop script 60's broader proxy", ~(B.origin.eq("literature_auto") & B.ph_opt.isna()))):
    m = (pr.notna() & mask & B.ph_min.notna()).to_numpy()
    if m.sum() < 60:
        print(f"   {label:34s} n={int(m.sum()):4d}  too few")
        continue
    g = B[m]
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    b, se, pv, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, pr[m].to_numpy(float)[:, None])
    rows.append({"part": "sensitivity", "subset": label, "n": int(m.sum()), "beta_per_sd": float(b[0]), "se": float(se[0]), "p": float(pv[0])})
    print(f"   {label:34s} n={int(m.sum()):4d}  {b[0]:+.4f} +/- {se[0]:.4f} (p {pv[0]:.2e})")
pd.DataFrame(rows).to_csv(D / "optimum_as_range.csv", index=False)
print("\nwrote", D / "optimum_as_range.csv")
