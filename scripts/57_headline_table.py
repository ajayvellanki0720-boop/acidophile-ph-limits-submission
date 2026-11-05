"""Recompute every headline figure in one place, with its units, panel, sample and model stated -> docs/NUMBERS.md and DATA_DIR/headline_numbers.csv
Five scale errors of the same family have been found in this project, and every one of them happened the same way: a number was quoted somewhere without its units or its
sample beside it, and later compared with a number computed on a different footing. The synthesis now quotes figures from two panels, three composites, several samples and
two scales, and a reader cannot tell which is which from the text alone. This script is the remedy and the regression test. It recomputes each figure from the raw tables in one
run, prints it on BOTH scales with its sample size and model, and checks the arithmetic that was got wrong before.
The check: for a predictor x, beta_per_unit * sd(x) must equal beta_per_sd. pgls_scan returns the per-SD figure, so the per-unit figure is derived by DIVIDING by sd(x), and
the script asserts the identity holds for every row rather than trusting it. Any future edit that reintroduces the confusion fails here.
Nothing in this script is a new result. If a number here disagrees with the synthesis, the synthesis is wrong.
Usage: python3 scripts/57_headline_table.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
WIDE_MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
WIDE_SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
HAND_MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
HAND_SOL = ["EF-Tu", "GroEL", "RecA"]
rl_w = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
rl_o = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size()
wide = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
orig = pd.read_csv(D / "structure_panel.csv").drop_duplicates("ncbi_assembly")
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))
SRCC = [f"src_{o}" for o in SRC]


def per_ref(tab, reflen, label, cls):
    c = f"{label}|{cls}|n"
    if c not in tab:
        return None
    m = tab.set_index("ncbi_assembly")
    n_ = m[c].reindex(B.ncbi_assembly)
    floor = max(8, 0.5 * reflen[(label, cls)]) if (label, cls) in reflen else 8
    num = (m[f"{label}|{cls}|D"] + m[f"{label}|{cls}|E"] - m[f"{label}|{cls}|K"] - m[f"{label}|{cls}|R"]).reindex(B.ncbi_assembly)
    return pd.Series((num / n_).where(n_ >= floor).to_numpy(), index=B.index)


def composite(tab, reflen, labels, cls, minrefs=3, cover=0.0, weights=None):
    cols = {}
    for r in labels:
        x = per_ref(tab, reflen, r, cls)
        if x is None or x.notna().mean() < cover:
            continue
        cols[r] = (x - x.mean()) / x.std()
    Zc = pd.DataFrame(cols)
    if weights == "pc":
        ccm = Zc.notna().all(axis=1)
        w = np.linalg.eigh(np.corrcoef(Zc[ccm].to_numpy(float), rowvar=False))[1][:, -1]
        if w.sum() < 0:
            w = -w
        W = pd.DataFrame(np.tile(w, (len(Zc), 1)), columns=Zc.columns, index=Zc.index).where(Zc.notna())
        s = (Zc * W).sum(axis=1, min_count=minrefs) / W.abs().sum(axis=1)
    else:
        s = Zc.mean(axis=1)
    return s.where(Zc.notna().sum(axis=1) >= minrefs), len(cols)


def fit(x, covs, require=()):
    """`require` restricts the sample to rows where those columns are present WITHOUT entering them as covariates, which is
    the only way to separate a change of model from a change of sample. The optimum is known for 1,655 of 2,768 genomes, so
    every comparison between an optimum-controlled row and a composition-only row otherwise mixes the two."""
    m = (x.notna() & B.ph_min.notna() & B[list(covs) + list(require)].notna().all(axis=1)).to_numpy()
    g = B[m]
    use = [c for c in covs + SRCC if g[c].nunique() > 1]
    X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
    v = x[m].to_numpy(float)
    b, se, p, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, v[:, None])
    sd = float(v.std())
    return {"n": int(m.sum()), "sd_of_predictor": sd, "beta_per_sd": float(b[0]), "se_per_sd": float(se[0]), "p": float(p[0]),
            "beta_per_unit": float(b[0]) / sd, "lam": float(lam), "covariates": ", ".join(["intercept"] + use),
            "sample_restricted_to": ", ".join(require) if require else ""}


SPECS = [
    ("exposed periplasmic, 32-reference rule-selected panel", wide, rl_w, WIDE_MEM, "peri_exp", 0.0, None, ["gw"], ()),
    # the pair below is the optimum decomposition on ONE sample, so the drop is a change of model alone
    ("exposed periplasmic, 32 refs, + composition, on the optimum-known sample", wide, rl_w, WIDE_MEM, "peri_exp", 0.0, None, ["gw"], ("ph_opt",)),
    ("exposed periplasmic, 32-reference panel, + growth optimum", wide, rl_w, WIDE_MEM, "peri_exp", 0.0, None, ["gw", "ph_opt"], ()),
    ("exposed periplasmic, 32 references, NO covariates", wide, rl_w, WIDE_MEM, "peri_exp", 0.0, None, [], ()),
    ("exposed periplasmic, 32 references, + optimum ONLY (no composition control)", wide, rl_w, WIDE_MEM, "peri_exp", 0.0, None, ["ph_opt"], ()),
    ("exposed periplasmic, 19 high-coverage references", wide, rl_w, WIDE_MEM, "peri_exp", 0.80, None, ["gw"], ()),
    ("exposed periplasmic, 19 high-coverage, factor-weighted", wide, rl_w, WIDE_MEM, "peri_exp", 0.80, "pc", ["gw"], ()),
    ("exposed periplasmic, 19 high-coverage, factor-weighted, + optimum", wide, rl_w, WIDE_MEM, "peri_exp", 0.80, "pc", ["gw", "ph_opt"], ()),
    ("exposed cytoplasmic control, 32 references", wide, rl_w, WIDE_MEM, "cyto_exp", 0.0, None, ["gw"], ()),
    ("transmembrane control, 32 references", wide, rl_w, WIDE_MEM, "TM", 0.0, None, ["gw"], ()),
    ("soluble-surface control, 5 references", wide, rl_w, WIDE_SOL, "surface", 0.0, None, ["gw"], ()),
    ("exposed periplasmic, 7-reference hand-picked panel", orig, rl_o, HAND_MEM, "peri_exp", 0.0, None, ["gw"], ()),
    ("exposed periplasmic, hand-picked panel, + growth optimum", orig, rl_o, HAND_MEM, "peri_exp", 0.0, None, ["gw", "ph_opt"], ()),
    ("soluble-surface control, hand-picked panel", orig, rl_o, HAND_SOL, "surface", 0.0, None, ["gw"], ()),
]
rows = []
for label, tab, reflen, labels, cls, cover, wmode, covs, require in SPECS:
    x, k = composite(tab, reflen, labels, cls, cover=cover, weights=wmode)
    if x is None or x.notna().sum() < 100:
        print(f"skipped (too few scorable): {label}")
        continue
    r = fit(x, covs, require)
    r.update({"figure": label, "references_used": k, "class": cls, "weighting": wmode or "unweighted mean"})
    rows.append(r)

rows_extra = []
R = pd.DataFrame(rows)
# the identity that five errors violated, asserted rather than assumed
bad = R[np.abs(R.beta_per_unit * R.sd_of_predictor - R.beta_per_sd) > 1e-9]
assert bad.empty, f"per-unit/per-SD identity failed for: {list(bad.figure)}"
print("per-unit x SD(predictor) == per-SD holds for all", len(R), "rows\n")
R.to_csv(D / "headline_numbers.csv", index=False)

print(f"{'figure':58s} {'refs':>4} {'n':>5} {'per SD':>9} {'per unit':>9} {'SD(x)':>6} {'p':>9}")
for x in R.itertuples():
    print(f"{x.figure[:58]:58s} {x.references_used:>4d} {x.n:>5d} {x.beta_per_sd:>+9.4f} {x.beta_per_unit:>+9.4f} {x.sd_of_predictor:>6.3f} {x.p:>9.1e}")

lines = ["# Headline numbers, with their units and samples",
         "",
         "Generated by `scripts/57_headline_table.py`. **If a figure in `SYNTHESIS.md` disagrees with this table, the synthesis is wrong.**",
         "",
         "Five scale errors of the same family have been found in this project, each one caused by a number quoted without its units or",
         "its sample beside it and later compared against a number computed on a different footing. Every figure below is given on both",
         "scales, with the panel, the composite, the covariates and the sample size that produced it. The script asserts",
         "`beta_per_unit * SD(predictor) == beta_per_sd` for every row, so an edit that reintroduces the confusion fails rather than",
         "publishing.",
         "",
         "`pgls_scan` returns the **per-standard-deviation** coefficient, because it standardises its predictor internally. The per-unit",
         "figure below is derived by DIVIDING by the predictor's SD. Multiplying the returned coefficient by the SD is the mistake.",
         "",
         "| Figure | Refs | n | per SD | per unit | SD(x) | p | lambda | Covariates |",
         "|---|---|---|---|---|---|---|---|---|"]
for x in R.itertuples():
    lines.append(f"| {x.figure} | {x.references_used} | {x.n:,} | **{x.beta_per_sd:+.4f}** | {x.beta_per_unit:+.4f} | {x.sd_of_predictor:.3f} | "
                 f"{x.p:.1e} | {x.lam:.3f} | {x.covariates} |")
lines += ["",
          "Notes on reading the table.",
          "",
          "- **Per SD is the figure to quote.** It is comparable across composites; per unit is not, because a composite averaging more",
          "  references has a smaller spread, so the same underlying effect gives a larger per-unit slope.",
          "- **The 32-reference and 7-reference panels are not the same measurement** and their per-unit figures should never be compared.",
          "- **The `+ growth optimum` rows use a smaller sample**, because the optimum is not known for every species. The drop from the",
          "  composition-only row to the optimum row therefore mixes a change of model with a change of sample.",
          "- The control rows exist so that the specificity claim can be checked on exactly the same footing as the claim itself.",
          ""]

# the optimum decomposition, on ONE sample, because the pair of rows above differ in model AND sample otherwise
try:
    same = R[R.figure.str.contains("on the optimum-known sample", regex=False)].iloc[0]
    with_opt = R[R.figure.str.contains("32-reference panel, + growth optimum", regex=False)].iloc[0]
    assert int(same.n) == int(with_opt.n), f"decomposition samples differ: {int(same.n)} vs {int(with_opt.n)}"
    share = 1 - with_opt.beta_per_sd / same.beta_per_sd
    print(f"\nThe growth-optimum decomposition, both rows on the same {int(same.n):,} genomes:")
    print(f"  composition only            {same.beta_per_sd:+.4f} per SD")
    print(f"  + growth optimum            {with_opt.beta_per_sd:+.4f} per SD")
    print(f"  so {share:.0%} of the association is attributable to the optimum, and {1 - share:.0%} is independent of it.")
    print(f"  Quoting the drop from the FULL-sample {float(R.iloc[0].beta_per_sd):+.4f} instead would mix a change of model with a")
    print(f"  change of sample, since the optimum is known for only {int(same.n):,} of {int(R.iloc[0].n):,} genomes.")
    lines += ["", "**The growth-optimum decomposition.** Both rows below are fitted on the same "
              f"{int(same.n):,} genomes, those whose optimum is known, so the drop is a change of model alone:",
              "",
              f"| Model | per SD |", "|---|---|",
              f"| composition only, optimum-known sample | {same.beta_per_sd:+.4f} |",
              f"| + growth optimum | {with_opt.beta_per_sd:+.4f} |",
              "",
              f"So **{share:.0%}** of the association is attributable to the growth optimum and **{1 - share:.0%}** is independent of it. "
              f"Comparing the optimum row against the full-sample {float(R.iloc[0].beta_per_sd):+.4f} would mix a change of model with a change of "
              f"sample and overstate the share, which is how the earlier \"roughly two-thirds\" figure arose.", ""]
except Exception as e:
    print("\n(decomposition not computed:", e, ")")


# DERIVED QUANTITIES. Correction 12 happened because a ratio quoted in prose was never a row in this table, so the table
# could not guard it. Anything this project states as a ratio, share or difference between two figures above is computed
# here instead of by hand, and a comparison whose two rows have unequal n is refused rather than printed, because that is
# the exact error correction 12 recorded.
def grab(fragment):
    """Exact match first, then a unique substring. An ambiguous fragment is an error, not a silent first-match -
    this assertion caught 'factor-weighted' matching both the plain and the optimum-controlled row on its first run."""
    exact = R[R.figure == fragment]
    if len(exact) == 1:
        return exact.iloc[0]
    hit = R[R.figure.str.contains(fragment, regex=False)]
    assert len(hit) == 1, f"{fragment!r} matched {len(hit)} rows: {list(hit.figure)}"
    return hit.iloc[0]


DERIVED = [
    ("share of the association attributable to the growth optimum", "32 refs, + composition, on the optimum-known sample",
     "32-reference panel, + growth optimum", "share", True),
    ("share attributable to genome-wide composition", "32 references, NO covariates",
     "32-reference rule-selected panel", "share", True),
    ("factor weighting against the unweighted mean, 19 high-coverage", "exposed periplasmic, 19 high-coverage, factor-weighted",
     "exposed periplasmic, 19 high-coverage references", "ratio", True),
    ("32-reference panel against the hand-picked 7", "32-reference rule-selected panel",
     "7-reference hand-picked panel", "ratio", False),
    ("exposed cytoplasmic control against the claim", "exposed cytoplasmic control",
     "32-reference rule-selected panel", "ratio", True),
    ("transmembrane control against the claim", "transmembrane control",
     "32-reference rule-selected panel", "ratio", True),
]
dlines = ["", "## Derived quantities", "",
          "Correction 12 happened because a ratio quoted in prose was never a row in this table, so the table could not guard it.",
          "Every ratio, share or difference this project states is computed here rather than by hand. A comparison whose two rows",
          "have unequal samples is refused rather than printed, since mixing a change of model with a change of sample is the",
          "error correction 12 recorded. **Per-unit figures are never compared across panels** — a composite averaging more",
          "references has a smaller spread, so the same underlying effect gives a larger per-unit slope.",
          "", "| Quantity | From | To | Value | Same sample |", "|---|---|---|---|---|"]
print("\nDerived quantities, computed rather than quoted:")
for name, a_frag, b_frag, kind, need_same in DERIVED:
    a, b = grab(a_frag), grab(b_frag)
    same = int(a.n) == int(b.n)
    if need_same and not same:
        print(f"  {name:58s} REFUSED: samples differ ({int(a.n):,} against {int(b.n):,})")
        dlines.append(f"| {name} | {a.figure} | {b.figure} | **refused**, samples differ ({int(a.n):,} vs {int(b.n):,}) | no |")
        rows_extra.append({"figure": name, "derived": True, "refused": True, "n_a": int(a.n), "n_b": int(b.n)})
        continue
    val = (1 - b.beta_per_sd / a.beta_per_sd) if kind == "share" else (a.beta_per_sd / b.beta_per_sd)
    txt = f"{val:.0%}" if kind == "share" else f"{val:+.2f}x"
    note = "yes" if same else "no, stated as a caution"
    print(f"  {name:58s} {txt:>8s}   (n {int(a.n):,} and {int(b.n):,})")
    dlines.append(f"| {name} | {a.figure} ({a.beta_per_sd:+.4f}) | {b.figure} ({b.beta_per_sd:+.4f}) | **{txt}** | {note} |")
    rows_extra.append({"figure": name, "derived": True, "kind": kind, "value": float(val), "n_a": int(a.n), "n_b": int(b.n),
                       "same_sample": same})
lines += dlines

(ROOT / "docs/NUMBERS.md").write_text("\n".join(lines))
pd.concat([R, pd.DataFrame(rows_extra)], ignore_index=True).to_csv(D / "headline_numbers.csv", index=False)
print("\nwrote", D / "headline_numbers.csv", "and", ROOT / "docs/NUMBERS.md")
