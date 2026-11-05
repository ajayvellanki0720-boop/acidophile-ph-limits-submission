"""Structure-based exposure over an EXPANDED reference panel -> DATA_DIR/structure_panel.csv (+ results/structure_panel_classes.csv).
Generalises 34a from three references to nine: membrane proteins YidC, SecY (Escherichia coli AND Bacillus subtilis, a reference-species control), ATP synthase subunit c,
NADH-quinone subunits L and N, cytochrome bo3 subunit 1; soluble controls EF-Tu, GroEL, RecA. ATP synthase subunit a was dropped: 80 % of its AlphaFold residues are below
pLDDT 70. Each position of a reference is classed from UniProt topology ('Extracellular' is treated as the periplasmic side) and from relative solvent accessibility of the
isolated model (Shrake-Rupley, Tien 2013 maxima; exposed = RSA >= 0.25). NEW relative to 34a: a position must have per-residue pLDDT >= 70 (PDB B-factor) to be classified at
all, since RSA of a poorly predicted residue is not interpretable; unclassified positions are reported as `lowconf` and excluded downstream. Every genome's best-aligning hit of
the Pfam family (BLOSUM62 local alignment) is tallied by the class of the reference position it aligns to; D, E, K, R and (since 2026-10-31) H are counted separately so
acidity, basicity and net charge can be separated later. Histidine is counted to test whether the residue ranking follows side-chain pKa rather than charge alone: see the
comment at the tally, and docs/SYNTHESIS.md section 8. Adding it is additive - no existing column changes - but the panel must be rebuilt for the new column to exist. Caveat: single-chain models, so lipid-facing transmembrane residues look exposed in isolation -> membrane proteins are classed by topology and
exposure is used only within the periplasmic and cytoplasmic segments.
Usage: python3 scripts/36a_structure_panel.py phenotypes.csv PROTEOME_DIR OUT.csv [processes]"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import pandas as pd
import pyhmmer
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
# label -> (uniprot, pfam, membrane?)
PANEL = {
    "YidC": ("P25714", "PF02096", True), "SecY": ("P0AGA2", "PF00344", True), "SecY_Bsub": ("P16336", "PF00344", True),
    "ATPsynC": ("P68699", "PF00137", True), "NuoL": ("P33607", "PF00361", True), "NuoN": ("P0AFF0", "PF00361", True), "CyoB": ("P0ABI8", "PF00115", True),
    "EF-Tu": ("P0CE47", "PF00009", False), "GroEL": ("P0A6F5", "PF00118", False), "RecA": ("P0A7G6", "PF00154", False),
}
# An optional wider panel (scripts/40a) can replace the built-in one; the default is unchanged so earlier results reproduce.
if os.environ.get("PANEL_CSV"):
    _r = pd.read_csv(ROOT / os.environ["PANEL_CSV"])
    _r = _r[_r.accepted]
    PANEL = {str(x.label): (str(x.uniprot), str(x.pfam), bool(x.membrane)) for x in _r.itertuples()}
MIN_PLDDT = 70.0
S = ROOT / "data/raw/structures"
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)
OUTSIDE = {"Periplasmic", "Extracellular"}


def _hmmpath(fam):
    for d in ("data/raw/hmm_markers", "data/raw/hmm", "data/raw/hmm_panel"):
        q = ROOT / d / f"{fam}.hmm"
        if q.exists():
            return q
    raise FileNotFoundError(fam)


def reference(acc, membrane):
    st = PDBParser(QUIET=True).get_structure(acc, S / f"AF-{acc}.pdb")
    ShrakeRupley().compute(st[0], level="R")
    seq, rsa, plddt = [], [], []
    for r in st[0]["A"]:
        aa = seq1(r.get_resname()); seq.append(aa); rsa.append(r.sasa / MAXASA[aa]); plddt.append(next(r.get_atoms()).get_bfactor())
    topo = ["?"] * len(seq)
    if membrane:
        for f in json.loads((S / f"{acc}.json").read_text())["features"]:
            if f["type"] not in ("Transmembrane", "Topological domain"):
                continue
            a, b = f["location"]["start"]["value"], f["location"]["end"]["value"]
            lab = "TM" if f["type"] == "Transmembrane" else ("peri" if f.get("description") in OUTSIDE else "cyto" if f.get("description") == "Cytoplasmic" else None)
            if lab:
                for i in range(a - 1, b): topo[i] = lab
    cls = []
    for t, x, p in zip(topo, rsa, plddt):
        if p < MIN_PLDDT:
            cls.append("lowconf")
        elif membrane:
            cls.append("TM" if t == "TM" else "unassigned" if t == "?" else f"{t}_{'exp' if x >= 0.25 else 'bur'}")
        else:
            cls.append("surface" if x >= 0.25 else "buried")
    return clean("".join(seq)), cls


REF = {}
HMMS = {}
ALN = Align.PairwiseAligner(); ALN.mode = "local"; ALN.substitution_matrix = substitution_matrices.load("BLOSUM62"); ALN.open_gap_score = -11; ALN.extend_gap_score = -1
ALPHA = set(ALN.substitution_matrix.alphabet)
clean = lambda s: "".join(c if c in ALPHA else "X" for c in s)   # selenocysteine (U) and other non-BLOSUM62 letters


def init():
    for name, (acc, fam, mem) in PANEL.items():
        p = _hmmpath(fam)
        with pyhmmer.plan7.HMMFile(p) as hf:
            HMMS[name] = next(iter(hf))
        REF[name] = reference(acc, mem)


def one(a):
    acc, d = a
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(Path(d) / f"{acc}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by = {_s(s.name): s for s in seqs}
    perfam = {}
    for name, hmm in HMMS.items():
        fam = PANEL[name][1]
        if fam not in perfam:
            perfam[fam] = {_s(h.name) for top in pyhmmer.hmmsearch([hmm], seqs, bit_cutoffs="gathering", cpus=1) for h in top if h.included}
    row = {"ncbi_assembly": acc}
    for name in PANEL:
        rseq, rcls = REF[name]
        best = None
        # DETERMINISM, correction 17. perfam[...] is a SET of protein names and Python randomises string
        # hashing per process, so this loop used to run in a different order in every run. With `>` below,
        # a TIE between two equally-aligning paralogues keeps whichever was seen first, which made the
        # output a function of the process rather than of the inputs: rebuilding the panel reproduced only
        # 99.71% of cells, and 31 of 32 references differed somewhere (scripts/79). Sorting makes the
        # tie-break the lexicographically first name, which is arbitrary but fixed.
        for nm in sorted(perfam[PANEL[name][1]]):
            q = clean(by[nm].textize().sequence.rstrip("*"))
            if len(q) <= 30:
                continue
            al = ALN.align(rseq, q)[0]
            if best is None or al.score > best[0]:
                best = (al.score, q, al)
        if best is None:
            continue
        score, q, al = best
        cnt = {}
        for (r0, r1), (q0, q1) in zip(*al.aligned):
            for i, j in zip(range(r0, r1), range(q0, q1)):
                c = rcls[i]
                if c in ("lowconf", "unassigned"):
                    continue
                aa = q[j]
                # HISTIDINE added 2026-10-31, purely additively: every existing column is unchanged and one new
                # `|H` column appears per class. It is here to test whether the residue ranking follows side-chain
                # pKa. Glutamate (pKa ~4.25) carries about 2.4x what aspartate (~3.65) does (scripts/75), which is
                # what a protonation account predicts - the carboxylate whose charge state turns over first as the
                # envelope acidifies is the one that matters. Histidine (~6.0) turns over earlier still, so a
                # charge-burden account predicts it behaves like a stronger glutamate. A BUFFERING account predicts
                # the opposite sign, since histidine is the main protein buffer, and this project's own withdrawn
                # histidine lead ran the other way proteome-wide. The prediction is therefore ambiguous, which is
                # exactly why the test discriminates between the two mechanisms instead of confirming either.
                v = cnt.setdefault(c, dict(n=0, D=0, E=0, K=0, R=0, H=0)); v["n"] += 1
                if aa in "DEKRH": v[aa] += 1
        row[f"{name}_score"] = score; row[f"{name}_len"] = len(q)
        for c, v in cnt.items():
            for k, x in v.items(): row[f"{name}|{c}|{k}"] = x
    return row


if __name__ == "__main__":
    phen, d, out = sys.argv[1:4]; procs = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    init()
    pd.DataFrame([{"ref": n, "pos": i + 1, "aa": REF[n][0][i], "class": REF[n][1][i]} for n in REF for i in range(len(REF[n][0]))]).to_csv(ROOT / os.environ.get("CLASSES_OUT", "results/structure_panel_classes.csv"), index=False)
    for n in REF: print(f"{n:10s}", pd.Series(REF[n][1]).value_counts().to_dict(), flush=True)
    accs = [a for a in dict.fromkeys(pd.read_csv(phen).ncbi_assembly.dropna()) if (Path(d) / f"{a}.faa").exists()]
    print(len(accs), "genomes", flush=True)
    with Pool(procs, initializer=init) as p:
        rows = list(p.imap_unordered(one, [(a, d) for a in accs], chunksize=8))
    # DETERMINISM, second source, found by the two-build verification of correction 17. The tie-break fix above made
    # every CELL reproducible, but the two builds still differed as FILES: imap_unordered yields a genome as whatever
    # worker finishes it, so row order - and with it the column order pandas infers from first appearance - varied per
    # run. Nothing downstream reads this file positionally (all 30 consumers merge on ncbi_assembly), so this changed
    # no result; the point is that "reproducible from its inputs" has to mean the file, not just the numbers in it.
    out_df = pd.DataFrame(rows).sort_values("ncbi_assembly").reset_index(drop=True)
    out_df = out_df[["ncbi_assembly"] + sorted(c for c in out_df.columns if c != "ncbi_assembly")]
    out_df.to_csv(out, index=False); print("done")
