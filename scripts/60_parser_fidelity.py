"""How good is the parser that produced 92% of the phenotypes, and does its one failure mode matter? -> DATA_DIR/parser_fidelity.csv
"The phenotype data are mostly parser-extracted" has sat in the synthesis's list of weaknesses since the dataset was broadened, never quantified. It is the weakness with the
most leverage, since 2,623 of 2,768 bacterial minima come from a rule-based extraction of species-description abstracts, so this puts a number on it. Three parts.
  1. FIDELITY against the hand-curated rows, from the validation the parser was built with (results/full/parser_validation.csv). This is reported here because it is a strong
     result that was never carried into the synthesis.
  2. THE GAP IN THAT VALIDATION. It was run on the abstracts that had been hand-curated, which were selected for being curatable and are acidophile-heavy: mean hand-read
     minimum 3.38 with a third of them at or below pH 3, against 5.94 and 0.5% in the parser rows actually in use. So it validates the parser on a different kind of abstract
     from most of the ones it was applied to, and a hand-check of 25 random rows from the range that carries the data (pH_min between 4 and 7) was done to close that gap.
     That check found 3 of 25 wrong, all with the same failure: a sentence that reports ONLY an optimum range, as in "grew optimally at 35 C, at pH 7.0-8.0", has its optimum
     recorded as the growth range. The error biases the minimum UPWARD, and in all three cases the signature is the same - no optimum was parsed, because the optimum had
     already been consumed as the range.
  3. WHETHER IT MATTERS. The signature gives a test. Rows with no parsed optimum (1,132 of 2,623) are the ones at risk; roughly a third of them by the hand-check, since many
     abstracts legitimately report a range and no optimum. If the failure mode mattered, dropping those rows should move the estimate. Both halves are fitted, and the
     hand-curated and database tiers are carried alongside as references that do not have this problem at all.
Usage: python3 scripts/60_parser_fidelity.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
rows = []

print("1. fidelity against the hand-curated rows")
v = pd.read_csv(ROOT / "results/full/parser_validation.csv")
v["err"] = v.ph_min_auto - v.ph_min_hand
for name, sub in (("all parsed rows", v), ("rows flagged clean", v[v.clean]), ("rows flagged usable (the gate in use)", v[v.usable_noopt])):
    s = sub.dropna(subset=["err"])
    rows.append({"part": "fidelity", "subset": name, "n": len(s), "min_exact_rate": float(sub.min_ok.mean()),
                 "mean_error": float(s.err.mean()), "sd_error": float(s.err.std()), "max_abs_error": float(s.err.abs().max())})
    print(f"   {name:38s} n={len(s):4d}  exact in {sub.min_ok.mean():.1%}  mean error {s.err.mean():+.3f}  SD {s.err.std():.3f}  "
          f"largest |error| {s.err.abs().max():.2f}")
print("   So extraction error is NOT the weakness: the parser reproduces what the abstract says almost exactly.")

print("\n2. the gap in that validation")
p = pd.read_csv(D / "phenotypes_all.csv")
B0 = p[p.domain == "Bacteria"].drop_duplicates("ncbi_assembly")
auto = B0[B0.origin == "literature_auto"]
print(f"   validation abstracts: mean hand-read minimum {v.ph_min_hand.mean():.2f}, {(v.ph_min_hand <= 3).mean():.0%} at or below pH 3")
print(f"   parser rows in use:   mean minimum {auto.ph_min.mean():.2f}, {(auto.ph_min <= 3).mean():.1%} at or below pH 3")
print("   hand-check of 25 random rows with pH_min in [4,7]: 3 wrong, all optimum-read-as-range, all with no parsed optimum")
rows.append({"part": "hand check", "subset": "25 random rows, pH_min in [4,7]", "n": 25, "errors": 3, "error_rate": 0.12,
             "note": "all three read an optimum range as the growth range; all three have no parsed optimum"})

print("\n3. does it matter? the failure mode's signature is a missing optimum")
has, no = auto.ph_opt.notna(), auto.ph_opt.isna()
print(f"   {int(has.sum())} parser rows have a parsed optimum, {int(no.sum())} do not")
print(f"   mean minimum {auto.ph_min[has].mean():.3f} with an optimum against {auto.ph_min[no].mean():.3f} without "
      f"(+{auto.ph_min[no].mean() - auto.ph_min[has].mean():.3f}, the direction the failure predicts)")
rows.append({"part": "signature", "subset": "parser rows with a parsed optimum", "n": int(has.sum()), "mean_ph_min": float(auto.ph_min[has].mean())})
rows.append({"part": "signature", "subset": "parser rows with no parsed optimum", "n": int(no.sum()), "mean_ph_min": float(auto.ph_min[no].mean())})

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

SUBSETS = [
    ("all bacteria", pd.Series(True, index=B.index)),
    ("drop parser rows with no parsed optimum", ~(B.origin.eq("literature_auto") & B.ph_opt.isna())),
    ("ONLY parser rows with no parsed optimum", B.origin.eq("literature_auto") & B.ph_opt.isna()),
    ("ONLY parser rows with a parsed optimum", B.origin.eq("literature_auto") & B.ph_opt.notna()),
    ("no parser rows at all (hand-curated + databases)", ~B.origin.eq("literature_auto")),
]
print(f"\n   exposed-periplasmic composite, genome-wide composition controlled (per SD)")
for label, mask in SUBSETS:
    m = (pr.notna() & mask & B.ph_min.notna()).to_numpy()
    if m.sum() < 60:
        print(f"   {label:50s} n={int(m.sum()):4d}  too few")
        continue
    g = B[m]
    use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    b, se, pv, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, pr[m].to_numpy(float)[:, None])
    rows.append({"part": "sensitivity", "subset": label, "n": int(m.sum()), "beta_per_sd": float(b[0]), "se": float(se[0]), "p": float(pv[0])})
    print(f"   {label:50s} n={int(m.sum()):4d}  {b[0]:+.4f} +/- {se[0]:.4f} (p {pv[0]:.2e})")

R = pd.DataFrame(rows)
R.to_csv(D / "parser_fidelity.csv", index=False)
print("\nwrote", D / "parser_fidelity.csv")
s = R[R.part == "sensitivity"].set_index("subset")
if "all bacteria" in s.index and "drop parser rows with no parsed optimum" in s.index:
    a, d_ = s.loc["all bacteria"], s.loc["drop parser rows with no parsed optimum"]
    print(f"\nDropping every row that could carry the failure mode moves the estimate from {a.beta_per_sd:+.4f} to {d_.beta_per_sd:+.4f} "
          f"({d_.beta_per_sd - a.beta_per_sd:+.4f}, {abs(d_.beta_per_sd - a.beta_per_sd) / a.se:.2f} standard errors of the full estimate).")
