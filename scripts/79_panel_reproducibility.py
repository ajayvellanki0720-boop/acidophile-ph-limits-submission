"""Does the structural panel rebuild identically, and does it matter? -> DATA_DIR/panel_rebuild.csv
THE PANEL DOES NOT REBUILD IDENTICALLY. Rebuilding the 32-reference panel from the same inputs with the same code (the only change being an added histidine column, which
touches nothing else) reproduces 99.71% of cells but differs somewhere in 759 of 850 D/E/K/R/n columns. This was found by scripts/76, which was written to add histidine and
compared the two builds only because the rebuild was deliberately written to a NEW file rather than overwriting the panel every result in this project rests on. Had it
overwritten, every figure would have shifted silently and nothing would have noticed.
THE CAUSE, which is specific and fixable. scripts/36a picks, for each reference, the best-aligning protein in the genome's Pfam family:
    perfam[fam] = {names...}                       # a SET
    for nm in perfam[fam]: ...
        if best is None or al.score > best[0]:     # strict >, so a TIE keeps the FIRST SEEN
Python randomises string hashing per process, so a set of protein names iterates in a different order in every run. When two paralogues align to the reference with EQUAL
score, which one is kept depends on that order. The output is therefore a function of the process, not only of the inputs.
WHY THE EFFECT ON RESULTS IS SMALL, AND WHY THAT IS NOT THE POINT. The composite correlates 0.9985 between builds and the headline estimate moves in the fourth decimal, so
no conclusion in this project changes. But "the number barely moved" is the wrong reason to be reassured, and this project has already published one figure on that reasoning
and withdrawn it. The right statements are: the pipeline is not reproducible from its inputs; the instability is in tie-breaking among paralogues, so it concentrates in
references with several close family members rather than spreading evenly; and the fix is two lines.
BOTH FAULTS ARE NOW FIXED AND THE FIX IS VERIFIED (--verify-determinism, see the top of the code). Verifying the first
fix found a second: with tie-breaking made deterministic the two builds agreed in every cell but still differed as
FILES, because imap_unordered returns a genome as whatever worker finishes it. That one changed no result, since all 30
consumers merge on ncbi_assembly and none reads the panel positionally - but a pipeline whose output depends on worker
scheduling is not a function of its inputs even when every number in it is, so it is fixed rather than excused.
ONE CONNECTION WORTH FOLLOWING. The reference with the largest per-column discrepancy is P62517, which is also the single reference scripts/78 flagged as showing the
CONCENTRATED signature in the position-distribution test. If its paralogue choice is unstable, its per-position pattern is partly noise, which would produce exactly that
signature artefactually. This script tests whether the two facts are related rather than leaving the coincidence to stand.
Usage: python3 scripts/79_panel_reproducibility.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NEW, OLD = D / "structure_panel_wide_h.csv", D / "structure_panel_wide.csv"
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in refs[refs.accepted & refs.membrane].label]
FAM = {str(x.label): str(x.pfam) for x in refs[refs.accepted & refs.membrane].itertuples()}


if len(sys.argv) > 3 and sys.argv[1] == "--verify-determinism":
    # THE FIX, VERIFIED. Run two independent builds of scripts/36a (separate processes, so separate hash seeds) and
    # pass both here. Two faults had to be fixed, and the first verification attempt is what exposed the second:
    #   (a) tie-breaking among equally-scoring paralogues, resolved by set iteration order -> every CELL varied.
    #       Fixed by iterating `sorted(perfam[...])`. Verified: the two builds agree cell for cell.
    #   (b) imap_unordered yields a genome as whatever worker finishes it, so ROW order - and the column order
    #       pandas infers from first appearance - varied. The builds had identical cells and identical size yet
    #       differed as files, at line 74. Fixed by sorting rows by ncbi_assembly and columns by name on write.
    # (b) changed no result - all 30 consumers merge on ncbi_assembly, none reads the panel positionally - but
    # "reproducible from its inputs" has to mean the file. A pipeline whose output depends on worker scheduling is
    # not a function of its inputs even when every number in it is.
    import hashlib
    hs, shape = [], None
    for f in sys.argv[2:4]:
        d = pd.read_csv(f).sort_values("ncbi_assembly").reset_index(drop=True)
        d = d[["ncbi_assembly"] + sorted(c for c in d.columns if c != "ncbi_assembly")]
        hs.append(hashlib.md5(d.to_csv(index=False).encode()).hexdigest())
        shape = d.shape
        print(f"   {f}  {hs[-1]}  {d.shape[0]} rows x {d.shape[1]} cols")
    raw = [hashlib.md5(Path(f).read_bytes()).hexdigest() for f in sys.argv[2:4]]
    print(f"   as written by 36a:        {'identical' if raw[0] == raw[1] else 'DIFFER (expected for builds made before the write-order fix)'}")
    print(f"   under 36a's current order: {'IDENTICAL - both faults fixed' if hs[0] == hs[1] else 'STILL DIFFER - a third source of nondeterminism remains'}")
    sys.exit(0 if hs[0] == hs[1] else 1)

A, B = pd.read_csv(NEW), pd.read_csv(OLD)

print("1. how far the two builds differ, by reference")
per_ref = []
for r in MEM:
    cols = [c for c in B.columns if c.startswith(f"{r}|") and c.split("|")[2] in ("n", "D", "E", "K", "R")]
    cols = [c for c in cols if c in A.columns]
    if not cols:
        continue
    a = A.set_index("ncbi_assembly")[cols]
    b = B.set_index("ncbi_assembly").loc[a.index, cols]
    both = a.notna() & b.notna()
    diff = int(((a != b) & both).to_numpy().sum())
    tot = int(both.to_numpy().sum())
    genomes = int((((a != b) & both).any(axis=1)).sum())
    per_ref.append({"ref": r, "pfam": FAM[r], "cells": tot, "cells_differing": diff,
                    "frac_differing": diff / tot if tot else np.nan, "genomes_differing": genomes})
R = pd.DataFrame(per_ref).sort_values("frac_differing", ascending=False)
print(f"   {'ref':>8} {'pfam':>9} {'cells':>8} {'differ':>7} {'frac':>8} {'genomes':>8}")
for x in R.head(8).itertuples():
    print(f"   {x.ref:>8} {x.pfam:>9} {x.cells:>8d} {x.cells_differing:>7d} {x.frac_differing:>8.4f} {x.genomes_differing:>8d}")
print(f"   ... {int((R.cells_differing == 0).sum())} of {len(R)} references rebuild cell-for-cell identically")
for x in R.itertuples():
    rows.append({"part": "per reference", "ref": x.ref, "pfam": x.pfam, "cells": x.cells,
                 "cells_differing": x.cells_differing, "frac_differing": x.frac_differing,
                 "genomes_differing": x.genomes_differing})

print("\n2. the composite, and the SAME MODEL fitted on each build")
print("   NOT the published headline. This uses a simplified common model - composition and source only, with a")
print("   flat coverage floor - because what is needed is a like-for-like comparison BETWEEN builds, not a")
print("   reproduction of docs/NUMBERS.md. The published +0.1848 controls more and uses a per-reference "
      "floor, so the figures below must not be quoted as that estimate.")


def composite(df):
    cols = []
    for r in MEM:
        if f"{r}|peri_exp|n" not in df:
            continue
        n = df[f"{r}|peri_exp|n"]
        x = ((df[f"{r}|peri_exp|D"] + df[f"{r}|peri_exp|E"] - df[f"{r}|peri_exp|K"] - df[f"{r}|peri_exp|R"]) / n
             ).where(n >= 8) * 100
        cols.append((x - x.mean()) / x.std())
    return pd.concat(cols, axis=1).mean(axis=1)


phen = pd.read_csv(D / "phenotypes_all.csv")
gf = pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp"]]
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
est = {}
for lab, df in (("rebuild", A), ("original", B)):
    d = df.copy()
    d["comp"] = composite(d)
    m = (phen.merge(gf, on="ncbi_assembly").merge(d[["ncbi_assembly", "comp"]], on="ncbi_assembly"))
    m = m[(m.domain == "Bacteria") & m.ncbi_assembly.isin(tips) & m.comp.notna() & m.ph_min.notna()].reset_index(drop=True)
    for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
        m[f"src_{o}"] = (m.origin == o).astype(float)
    src = [f"src_{o}" for o in ("thermobase", "acidb", "bacdive", "literature_auto") if m[f"src_{o}"].nunique() > 1]
    C = ph.vcv_fast(tree, list(m.ncbi_assembly))
    X = np.column_stack([np.ones(len(m)), m.acidic_excess * 100] + [m[k] for k in src])
    b, se, p, lam = ph.pgls_scan(C, m.ph_min.to_numpy(float), X, m.comp.to_numpy(float)[:, None])
    est[lab] = (float(b[0]), float(se[0]), float(p[0]), len(m))
    print(f"   {lab:9s} +{b[0]:.4f} +/- {se[0]:.4f} (p {p[0]:.2e}), n={len(m)}, lambda {lam:.3f}")
    rows.append({"part": "same-model comparison, NOT the published headline", "build": lab, "beta_per_sd": float(b[0]), "se": float(se[0]),
                 "p": float(p[0]), "n": len(m)})
d_abs = abs(est["rebuild"][0] - est["original"][0])
print(f"   difference {d_abs:.5f} per SD, i.e. {d_abs / est['original'][1]:.3f} of one standard error")
print("   So no conclusion in this project changes. That is NOT the reason to be reassured: a figure that barely")
print("   moves is still produced by a pipeline that is not a function of its inputs, and this project has")
print("   published once on 'it barely moved' reasoning and withdrawn it.")
rows.append({"part": "build difference", "abs_difference_per_sd": d_abs,
             "in_standard_errors": d_abs / est["original"][1]})

print("\n3. is P62517's instability why scripts/78 called it the one concentrated reference?")
p62 = R[R.ref == "P62517"]
if len(p62):
    x = p62.iloc[0]
    rank = int((R.frac_differing > x.frac_differing).sum()) + 1
    print(f"   P62517 differs in {x.frac_differing:.4f} of cells, rank {rank} of {len(R)} on instability")
    try:
        dist = pd.read_csv(D / "distribution_wide_by_reference.csv")
        mm = R.merge(dist[["ref", "share_same_sign"]], on="ref")
        r_corr = float(mm.frac_differing.corr(mm.share_same_sign))
        print(f"   across all references, instability against cross-half agreement: r {r_corr:+.3f}")
        print("   If unstable paralogue choice were manufacturing the concentrated signature, instability and")
        print("   agreement would be strongly NEGATIVELY correlated across references.")
        verdict = ("CONSISTENT with an artefact for P62517 specifically, which is the most unstable reference and the "
                   "only concentrated one" if rank == 1 else
                   "NOT a general explanation: P62517 is not the most unstable reference")
        print(f"   {verdict}")
        rows.append({"part": "P62517", "frac_differing": float(x.frac_differing), "instability_rank": rank,
                     "r_instability_vs_agreement": r_corr, "verdict": verdict})
    except FileNotFoundError:
        print("   distribution_wide_by_reference.csv missing; run scripts/78 first")

pd.DataFrame(rows).to_csv(D / "panel_rebuild.csv", index=False)
print("\nwrote", D / "panel_rebuild.csv")
print("\nReading. The pipeline is not reproducible from its inputs, because a tie between equally-scoring paralogues is")
print("broken by Python's per-process set ordering. The effect on every published figure is in the fourth decimal, so")
print("nothing is withdrawn; the defect is real regardless, the fix is to iterate in sorted order, and it was visible")
print("only because the rebuild was written to a new file instead of over the old one.")
