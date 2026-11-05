"""What exactly tracks pH_min on exposed periplasmic surfaces, and is it distributed? -> DATA_DIR/exposure_decomposition.csv, DATA_DIR/exposure_decomposition_pooled.csv, DATA_DIR/exposure_position_splits.csv, DATA_DIR/exposure_position_splits_summary.csv, and with the `wide` argument DATA_DIR/exposure_decomposition_wide.csv, DATA_DIR/exposure_decomposition_pooled_wide.csv, DATA_DIR/exposure_position_splits_wide.csv, DATA_DIR/exposure_position_splits_summary_wide.csv.
Every output is declared here because scripts/69 builds a claim's candidate pool from this line: the pooled file was added on 2026-10-31 and left undeclared, so the two
figures docs/PREREGISTRATION.md quotes from it had no pool and the check reported them untraceable. Both panels' spellings are listed because that script keeps a panel
variant's pool separate from the default run's, which is what stops one panel's figures passing as the other's.
Part 1 (composition): for the exposed periplasmic class of each membrane reference of the 36a panel, PGLS of pH_min on each of D, E, D+E, K, R, K+R and the net excess
(D+E-K-R), all per percentage point and all with genome-wide acidic excess as a covariate, then random-effects pooling across references. This separates acidity from
basicity and aspartate from glutamate.
Part 2 (distribution): the reference positions of a class are split at random into two halves, acidic excess is recomputed from each half alone, and both halves are tested.
If the signal is distributed over many positions both halves should show it; if a few positions carry it, one half will and the other will not. Reported over SPLITS random
splits as the share of splits in which both halves are positive and the median beta of each half. Random halves are used rather than dropping the strongest positions,
because ranking positions by their own effect on the same data and re-testing would be circular.
PANEL ARGUMENT, added 2026-10-31. Both parts had only ever been run on the 7-reference hand-picked panel while the headline rests on the 32-reference
rule-selected one, so the two claims this script supplies to docs/SYNTHESIS.md - which residues carry the signal, and whether it is distributed over many
positions - were the last rows of that document's evidence table still quoted for the hand-picked panel alone. They are also the two the charge-reversal
experiment leans on most directly: it reverses D/E to K at fourteen positions at once, which assumes both that the acidic side is what matters and that no
small set of positions carries the effect. `wide` runs the same two tests on the rule-selected panel. Outputs take a _wide suffix so the hand-picked results
stay on disk for comparison rather than being overwritten.
Usage: python3 scripts/36c_exposure_decomposition.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import pyhmmer
from multiprocessing import Pool
from Bio import Align, Phylo
from Bio.Align import substitution_matrices
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
WHICH = (sys.argv[1] if len(sys.argv) > 1 else "orig").lower()
assert WHICH in ("orig", "wide"), f"panel must be orig or wide, got {WHICH!r}"
SPLITS = 40
if WHICH == "wide":
    TAG, PANEL_CSV, CLS_CSV = "_wide", "structure_panel_wide.csv", "results/structure_panel_wide_classes.csv"
    _r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
    _r = _r[_r.accepted & _r.membrane]
    MEM = [str(x) for x in _r.label]
    # Pfam per reference, needed by part 2 to find the orthologue in each proteome. The hand-picked panel's map was
    # written out by hand because it had four entries; this one is read from the same table the panel is defined by.
    PFAM = {str(x.label): str(x.pfam) for x in _r.itertuples()}
else:
    TAG, PANEL_CSV, CLS_CSV = "", "structure_panel.csv", "results/structure_panel_classes.csv"
    MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
    PFAM = {"YidC": "PF02096", "NuoL": "PF00361", "NuoN": "PF00361", "CyoB": "PF00115"}
print(f"panel: {WHICH} ({len(MEM)} membrane references), outputs tagged {TAG!r}")
CLS = pd.read_csv(ROOT / CLS_CSV)
reflen = CLS.groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / PANEL_CSV), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]


def fit(mask, xv, extra=("gw",)):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    # Source indicators are screened for constancy on the FULL sample in SRCC, but a fit runs on the subset of genomes
    # scorable for one reference, and a narrow reference can leave an indicator all-zero there - a singular design.
    # The hand-picked panel never hit this because all 7 of its references score nearly every genome; the 32-reference
    # panel does. The screen is the same one, applied to the sample actually being fitted.
    keep = [k for k in SRCC if g[k].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[k] for k in list(extra) + keep])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, np.asarray(xv, float)[:, None])
    return b[0], se[0], p[0], len(g)


def pool(bs, ses):
    bs, ses = np.asarray(bs), np.asarray(ses)
    w = 1 / ses ** 2; fe = (w * bs).sum() / w.sum(); Q = (w * (bs - fe) ** 2).sum(); k = len(bs) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (ses ** 2 + tau2); mu = (wr * bs).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem))


# ---------------- Part 1: composition (inside part1(), so spawned workers do not re-run it on import) ----------------
def part1():
    COMPS = {"D": ("D",), "E": ("E",), "D+E": ("D", "E"), "K": ("K",), "R": ("R",), "K+R": ("K", "R"), "net D+E-K-R": ("D", "E", "-K", "-R")}
    rows = []
    for ref in MEM:
        c = "peri_exp"
        if f"{ref}|{c}|n" not in B:
            continue
        n = B[f"{ref}|{c}|n"]; need = max(8, 0.5 * reflen[(ref, c)])
        for lab, parts in COMPS.items():
            num = sum((-1 if p.startswith("-") else 1) * B[f"{ref}|{c}|{p.lstrip('-')}"] for p in parts)
            x = (num / n).where(n >= need) * 100
            m = x.notna().to_numpy()
            if m.sum() < 200:
                continue
            b, se, p_, nn = fit(m, x[m].to_numpy())
            rows.append({"ref": ref, "component": lab, "n": nn, "sd_pts": x[m].std(), "beta_per_pt": b, "se": se, "p": p_})
    R1 = pd.DataFrame(rows); R1.to_csv(D / f"exposure_decomposition{TAG}.csv", index=False)
    pd.set_option("display.width", 200)
    print("Composition of the exposed periplasmic signal (per percentage point, genome-wide acidic excess controlled)")
    print(R1.round(4).to_string(index=False))
    print("\n  pooled across the membrane references (random effects):")
    # The pooled figures are what docs/SYNTHESIS.md quotes, so they are STORED as rows rather than only printed.
    # A derived quantity quoted from stdout is unguarded - correction 12 - and scripts/69 cannot trace it to a file.
    pooled = []
    for lab in COMPS:
        y = R1[R1.component == lab]
        if len(y) >= 2:
            mu, sem, pv = pool(y.beta_per_pt, y.se)
            print(f"    {lab:12s} {mu:+.3f} +/- {sem:.3f} (p {pv:.1e}); positive in {(y.beta_per_pt > 0).sum()}/{len(y)}")
            pooled.append({"panel": WHICH, "component": lab, "pooled_beta_per_pt": float(mu), "pooled_se": float(sem),
                           "pooled_p": float(pv), "n_refs": int(len(y)), "n_refs_positive": int((y.beta_per_pt > 0).sum()),
                           "median_beta_per_pt": float(y.beta_per_pt.median())})
    pd.DataFrame(pooled).to_csv(D / f"exposure_decomposition_pooled{TAG}.csv", index=False)
    print(f"  wrote {D / f'exposure_decomposition_pooled{TAG}.csv'}")

# ---------------- Part 2: distribution over positions ----------------
# Residue aligned to each exposed-periplasmic reference position, per genome (each proteome read once, genomes in parallel).
if WHICH == "wide":
    # On the wide panel the reference label IS its accession, so part 2 covers every membrane reference with enough
    # exposed periplasmic positions to split in half at all. A reference with fewer than 8 cannot support the test:
    # each half would be scored from under 4 aligned positions, which is the per-half floor used below.
    _npos = CLS[CLS["class"] == "peri_exp"].groupby("ref").size()
    MAXS = {r: (r, PFAM[r]) for r in MEM if int(_npos.get(r, 0)) >= 8}
else:
    MAXS = {"YidC": ("P25714", "PF02096"), "NuoL": ("P33607", "PF00361"), "NuoN": ("P0AFF0", "PF00361"), "CyoB": ("P0ABI8", "PF00115")}
ALN = Align.PairwiseAligner(); ALN.mode = "local"; ALN.substitution_matrix = substitution_matrices.load("BLOSUM62"); ALN.open_gap_score = -11; ALN.extend_gap_score = -1
ALPHA = set(ALN.substitution_matrix.alphabet)
clean = lambda s: "".join(ch if ch in ALPHA else "X" for ch in s)
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)
RS, POS, HMMS = {}, {}, {}


def init2():
    for ref, (acc, fam) in MAXS.items():
        RS[ref] = clean("".join(CLS[CLS.ref == ref].sort_values("pos").aa))
        POS[ref] = CLS[(CLS.ref == ref) & (CLS["class"] == "peri_exp")].pos.to_numpy() - 1
        q = next(z for d in ("data/raw/hmm_markers", "data/raw/hmm", "data/raw/hmm_panel") if (z := ROOT / d / f"{fam}.hmm").exists())
        with pyhmmer.plan7.HMMFile(q) as hf:
            HMMS[ref] = next(iter(hf))


def residues(a):
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(ROOT / "data/raw/proteomes" / f"{a}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by = {_s(s.name): s for s in seqs}
    out = {}
    hits = {}
    for ref in MAXS:
        fam = MAXS[ref][1]
        if fam not in hits:
            hits[fam] = {_s(h.name) for top in pyhmmer.hmmsearch([HMMS[ref]], seqs, bit_cutoffs="gathering", cpus=1) for h in top if h.included}
        best = None
        # DETERMINISM, correction 17, extended. The same unsorted-set tie-break as scripts/36a: with strict `>`
        # below, two equally-aligning paralogues resolve by Python's per-process set order. Sorting makes the
        # tie-break the lexicographically first name - arbitrary, but a function of the inputs.
        for nm in sorted(hits[fam]):
            q = clean(by[nm].textize().sequence.rstrip("*"))
            if len(q) <= 30:
                continue
            al = ALN.align(RS[ref], q)[0]
            if best is None or al.score > best[0]:
                best = (al.score, q, al)
        if best is None:
            continue
        score, q, al = best
        v = np.full(len(RS[ref]), "-", dtype="<U1")
        for (r0, r1), (q0, q1) in zip(*al.aligned):
            for i, j in zip(range(r0, r1), range(q0, q1)):
                v[i] = q[j]
        out[ref] = v[POS[ref]]
    return a, out


if __name__ == "__main__":
    part1()
    init2()
    with Pool(8, initializer=init2) as pl:
        got = dict(pl.imap_unordered(residues, list(B.ncbi_assembly), chunksize=8))
    rows2, summ2 = [], []
    for ref in MAXS:
        pos = POS[ref]
        arr = np.array([got.get(a, {}).get(ref, np.full(len(pos), "-", dtype="<U1")) for a in B.ncbi_assembly])
        rng = np.random.default_rng(20261025)
        for s in range(SPLITS):
            perm = rng.permutation(len(pos)); halves = (perm[: len(pos) // 2], perm[len(pos) // 2:])
            out = []
            for h in halves:
                sl = arr[:, h]
                nn = (sl != "-").sum(1)
                num = np.isin(sl, ["D", "E"]).sum(1) - np.isin(sl, ["K", "R"]).sum(1)
                x = pd.Series(np.where(nn >= max(4, 0.5 * len(h)), 100 * num / np.maximum(nn, 1), np.nan))
                m = x.notna().to_numpy()
                out.append(fit(m, x[m].to_numpy()) if m.sum() >= 200 else (np.nan,) * 4)
            rows2.append({"ref": ref, "split": s, "positions": len(pos), "beta_a": out[0][0], "p_a": out[0][2], "beta_b": out[1][0], "p_b": out[1][2]})
        z = pd.DataFrame(rows2); z = z[z.ref == ref]
        print(f"\n{ref} exposed periplasmic, {len(pos)} positions, {SPLITS} random half-splits: both halves positive in {((z.beta_a > 0) & (z.beta_b > 0)).mean():.0%} of splits | "
              f"median beta {np.nanmedian(z.beta_a):+.3f} and {np.nanmedian(z.beta_b):+.3f} | both p<0.05 in {((z.p_a < 0.05) & (z.p_b < 0.05)).mean():.0%}")
        # the SHARE is the quoted quantity and is derived from the per-split rows, so it is stored in its own column
        summ2.append({"panel": WHICH, "ref": ref, "positions": int(len(pos)), "n_splits": int(SPLITS),
                      "share_both_halves_positive": float(((z.beta_a > 0) & (z.beta_b > 0)).mean()),
                      "share_both_halves_p05": float(((z.p_a < 0.05) & (z.p_b < 0.05)).mean()),
                      "median_beta_half_a": float(np.nanmedian(z.beta_a)), "median_beta_half_b": float(np.nanmedian(z.beta_b)),
                      "n_splits_evaluable": int((z.beta_a.notna() & z.beta_b.notna()).sum())})
    pd.DataFrame(rows2).to_csv(D / f"exposure_position_splits{TAG}.csv", index=False)
    S2 = pd.DataFrame(summ2)
    S2.to_csv(D / f"exposure_position_splits_summary{TAG}.csv", index=False)
    if len(S2) > 2:
        print(f"\nacross {len(S2)} references: median share with both halves positive "
              f"{S2.share_both_halves_positive.median():.0%}, range {S2.share_both_halves_positive.min():.0%} to "
              f"{S2.share_both_halves_positive.max():.0%}; "
              f"{int((S2.share_both_halves_positive >= 0.8).sum())} of {len(S2)} at 80% or above")
    print("\nwrote", D / f"exposure_position_splits{TAG}.csv", "and", D / f"exposure_position_splits_summary{TAG}.csv")
