"""Per-position analysis of the exposed periplasmic classes -> DATA_DIR/candidate_positions.csv (+ DATA_DIR/position_residues_<ref>.csv.gz).
For each exposed periplasmic reference position of YidC, NuoL and NuoN: the reference (Escherichia coli) residue, how often the aligned residue is acidic or basic across bacteria,
the mean pH_min of the genomes that carry an acidic versus a basic residue there, and a PGLS of pH_min on an acidic-residue indicator for that single position (genome-wide
acidic excess controlled, source indicators included). Positions are then scored as mutagenesis targets: the reference must carry D or E (so a reversal is possible), the site
must be naturally polymorphic between acidic and basic across species (so both residues are known to fold and function somewhere), and the acidic state must be associated with
a higher pH_min. Naturally polymorphic sites are preferred over strictly conserved ones precisely because a conserved acidic residue is more likely to be structurally or
catalytically required, and mutating it would test folding rather than the hypothesis.
Usage: python3 scripts/37a_candidate_positions.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import json
import os
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import pyhmmer
from Bio import Align, Phylo
from Bio.Align import substitution_matrices
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
REFS = {"YidC": ("P25714", "PF02096"), "NuoL": ("P33607", "PF00361"), "NuoN": ("P0AFF0", "PF00361")}
if os.environ.get("REFS"):           # e.g. REFS='{"P69681": ["P69681", "PF00909"]}'
    REFS = {k: tuple(v) for k, v in json.loads(os.environ["REFS"]).items()}
# REFS and the class table are configurable so a reference can be added without duplicating this scan. The
# pre-registration recommends AmtB as a second forward arm, and scripts/71 runs exactly this code for it.
CLS = pd.read_csv(ROOT / os.environ.get("CLASS_CSV", "results/structure_panel_classes.csv"))
OUT_TAG = os.environ.get("OUT_TAG", "")
ALN = Align.PairwiseAligner(); ALN.mode = "local"; ALN.substitution_matrix = substitution_matrices.load("BLOSUM62"); ALN.open_gap_score = -11; ALN.extend_gap_score = -1
ALPHA = set(ALN.substitution_matrix.alphabet)
clean = lambda s: "".join(ch if ch in ALPHA else "X" for ch in s)
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)
RS, HMMS = {}, {}


def init():
    for ref, (acc, fam) in REFS.items():
        RS[ref] = clean("".join(CLS[CLS.ref == ref].sort_values("pos").aa))
        q = next(z for d in ("data/raw/hmm_markers", "data/raw/hmm", "data/raw/hmm_panel") if (z := ROOT / d / f"{fam}.hmm").exists())
        with pyhmmer.plan7.HMMFile(q) as hf:
            HMMS[ref] = next(iter(hf))


def residues(a):
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(ROOT / "data/raw/proteomes" / f"{a}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by = {_s(s.name): s for s in seqs}
    out, hits = {}, {}
    for ref, (acc, fam) in REFS.items():
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
        out[ref] = v
    return a, out


if __name__ == "__main__":
    df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly"))
    tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
    B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
    SRC = [f"src_{o}" for o in ("thermobase", "acidb", "bacdive", "literature_auto")]
    for c in SRC:
        B[c] = (B.origin == c[4:]).astype(float)
    Call = ph.vcv_fast(tree, list(B.ncbi_assembly)); B["gw"] = B.acidic_excess * 100
    SRCC = [c for c in SRC if B[c].nunique() > 1]
    init()
    with Pool(8, initializer=init) as pl:
        got = dict(pl.imap_unordered(residues, list(B.ncbi_assembly), chunksize=8))
    rows = []
    for ref in REFS:
        full = np.array([got.get(a, {}).get(ref, np.full(len(RS[ref]), "-", dtype="<U1")) for a in B.ncbi_assembly])
        pd.DataFrame(full, index=B.ncbi_assembly).to_csv(D / f"position_residues_{ref}.csv.gz", index=True)
        sub = CLS[(CLS.ref == ref) & (CLS["class"] == "peri_exp")]
        for _, r in sub.iterrows():
            col = full[:, int(r.pos) - 1]
            acid, base = np.isin(col, ["D", "E"]), np.isin(col, ["K", "R"])
            m = (acid | base)
            if m.sum() < 150 or acid.sum() < 25 or base.sum() < 25:
                keep = dict(beta=np.nan, se=np.nan, p=np.nan)
            else:
                g = B[m].reset_index(drop=True); Cm = Call[np.ix_(m, m)]
                X = np.column_stack([np.ones(len(g)), g.gw] + [g[c] for c in SRCC if g[c].nunique() > 1])
                b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, acid[m].astype(float)[:, None])
                keep = dict(beta=b[0], se=se[0], p=p[0])
            rows.append({"ref": ref, "pos": int(r.pos), "ref_aa": r.aa, "n_aligned": int((col != "-").sum()),
                         "frac_acidic": float(acid.mean()), "frac_basic": float(base.mean()),
                         "ph_min_acidic": float(B.ph_min[acid].mean()) if acid.sum() else np.nan,
                         "ph_min_basic": float(B.ph_min[base].mean()) if base.sum() else np.nan, **keep})
    R = pd.DataFrame(rows)
    R["polymorphic"] = (R.frac_acidic >= 0.15) & (R.frac_basic >= 0.15)
    R["target"] = R.ref_aa.isin(["D", "E"]) & R.polymorphic & (R.beta > 0)
    R = R.sort_values(["target", "beta"], ascending=[False, False])
    R.to_csv(D / f"candidate_positions{OUT_TAG}.csv", index=False)
    pd.set_option("display.width", 220)
    print(f"{R.target.sum()} candidate targets of {len(R)} exposed periplasmic positions "
          f"(reference carries D/E, acidic/basic polymorphic at >=15 % each, acidic state associated with higher pH_min)\n")
    print(R[R.target].head(25)[["ref", "pos", "ref_aa", "n_aligned", "frac_acidic", "frac_basic", "ph_min_acidic", "ph_min_basic", "beta", "se", "p"]].round(3).to_string(index=False))
    print("\nper reference: candidates /", R.groupby("ref").target.agg(["sum", "size"]).to_dict("index"))
