"""Structure-based exposure of acidic residues in conserved proteins -> DATA_DIR/structure_exposure.csv (+ results/structure_reference_classes.csv).
Reference: AlphaFold v6 models of E. coli YidC (P25714), SecY (P0AGA2), EF-Tu (P0CE47; soluble control), plus UniProt topology. Each residue of the reference is classed from
UniProt topology (membrane proteins: periplasmic / cytoplasmic / transmembrane) and from relative solvent accessibility of the isolated model (Shrake-Rupley, Tien 2013 maxima;
exposed = RSA >= 0.25). Every genome's best-matching hit of the Pfam family is aligned locally to the reference (BLOSUM62) and its residues are tallied by the class of the
reference position they align to. Caveat: models are single chains; lipid-facing TM residues look 'exposed' in isolation, so membrane proteins are classed by topology,
with exposure used only inside the periplasmic and cytoplasmic segments.
Usage: python3 scripts/34a_structure_exposure.py phenotypes.csv PROTEOME_DIR OUT.csv [processes]"""
import json
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
REFS = {"YidC": ("P25714", "PF02096", True), "SecY": ("P0AGA2", "PF00344", True), "EF-Tu": ("P0CE47", "PF00009", False)}
S = ROOT / "data/raw/structures"
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)


def reference(acc, membrane):
    st = PDBParser(QUIET=True).get_structure(acc, S / f"AF-{acc}.pdb")
    ShrakeRupley().compute(st[0], level="R")
    seq, rsa = [], []
    for r in st[0]["A"]:
        aa = seq1(r.get_resname()); seq.append(aa); rsa.append(r.sasa / MAXASA[aa])
    topo = ["?"] * len(seq)
    if membrane:
        for f in json.loads((S / f"{acc}.json").read_text())["features"]:
            a, b = f["location"]["start"]["value"], f["location"]["end"]["value"]
            lab = {"Transmembrane": "TM"}.get(f["type"]) or {"Periplasmic": "peri", "Cytoplasmic": "cyto"}.get(f.get("description", ""))
            if f["type"] in ("Transmembrane", "Topological domain") and lab:
                for i in range(a - 1, b): topo[i] = lab
    cls = []
    for t, x in zip(topo, rsa):
        if membrane:
            cls.append("TM" if t == "TM" else f"{t}_{'exp' if x >= 0.25 else 'bur'}")
        else:
            cls.append("surface" if x >= 0.25 else "buried")
    return "".join(seq), cls


REF = {}
ALN = Align.PairwiseAligner(); ALN.mode = "local"; ALN.substitution_matrix = substitution_matrices.load("BLOSUM62"); ALN.open_gap_score = -11; ALN.extend_gap_score = -1
HMMS = {}


def init():
    for name, (_, fam, _) in REFS.items():
        with pyhmmer.plan7.HMMFile(ROOT / "data/raw/hmm_markers" / f"{fam}.hmm") as hf:
            HMMS[name] = next(iter(hf))
    for name, (acc, _, mem) in REFS.items():
        REF[name] = reference(acc, mem)


def one(a):
    acc, d = a
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(Path(d) / f"{acc}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by = {_s(s.name): s for s in seqs}
    row = {"ncbi_assembly": acc}
    for name, hmm in HMMS.items():
        names = {_s(h.name) for top in pyhmmer.hmmsearch([hmm], seqs, bit_cutoffs="gathering", cpus=1) for h in top if h.included}
        rseq, rcls = REF[name]
        best = None
        # DETERMINISM, correction 17, extended. The same unsorted-set tie-break as scripts/36a: with strict `>`
        # below, two equally-aligning paralogues resolve by Python's per-process set order. Sorting makes the
        # tie-break the lexicographically first name - arbitrary, but a function of the inputs.
        for nm in sorted(names):
            q = by[nm].textize().sequence.rstrip("*")
            al = ALN.align(rseq, q)[0] if len(q) > 30 else None
            if al is not None and (best is None or al.score > best[0]): best = (al.score, q, al)
        if best is None:
            continue
        score, q, al = best
        cnt = {}
        for (r0, r1), (q0, q1) in zip(*al.aligned):
            for i, j in zip(range(r0, r1), range(q0, q1)):
                c = rcls[i]; aa = q[j]
                d_ = cnt.setdefault(c, dict(n=0, D=0, E=0, K=0, R=0)); d_["n"] += 1
                if aa in "DEKR": d_[aa] += 1
        row[f"{name}_score"] = score; row[f"{name}_len"] = len(q)
        for c, v in cnt.items():
            for k, x in v.items(): row[f"{name}|{c}|{k}"] = x
    return row


if __name__ == "__main__":
    phen, d, out = sys.argv[1:4]; procs = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    init()
    pd.DataFrame([{"ref": n, "pos": i + 1, "aa": REF[n][0][i], "class": REF[n][1][i]} for n in REF for i in range(len(REF[n][0]))]).to_csv(ROOT / "results/structure_reference_classes.csv", index=False)
    for n in REF: print(n, pd.Series(REF[n][1]).value_counts().to_dict(), flush=True)
    accs = [a for a in dict.fromkeys(pd.read_csv(phen).ncbi_assembly.dropna()) if (Path(d) / f"{a}.faa").exists()]
    print(len(accs), "genomes", flush=True)
    with Pool(procs, initializer=init) as p:
        rows = list(p.imap_unordered(one, [(a, d) for a in accs], chunksize=8))
    # DETERMINISM, correction 17, extended. imap_unordered yields a genome as whatever worker finishes it, so row
    # order - and the column order pandas infers from first appearance - varied per run. 34a is the three-reference
    # ancestor of 36a and carried the same two faults; nothing reads it positionally, so no figure changes.
    out_df = pd.DataFrame(rows).sort_values("ncbi_assembly").reset_index(drop=True)
    out_df = out_df[["ncbi_assembly"] + sorted(c for c in out_df.columns if c != "ncbi_assembly")]
    out_df.to_csv(out, index=False); print("done")
