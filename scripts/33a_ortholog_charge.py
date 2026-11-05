"""Charge composition of ORTHOLOGOUS conserved proteins (same Pfam family in every genome), membrane vs soluble -> DATA_DIR/ortholog_charge.csv.
Residues of every gathering-cutoff hit of a family are pooled per genome. Usage: python3 scripts/33a_ortholog_charge.py phenotypes.csv PROTEOME_DIR OUT.csv [processes]"""
import sys
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
import pyhmmer
from acidophile.proteome import AA

FAMILIES = {"membrane": ["PF00344", "PF00119", "PF00137", "PF00115", "PF00361", "PF02096"], "soluble": ["PF00009", "PF00154", "PF00118", "PF00006", "PF00012"]}
HMMS = None
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)


def init():
    global HMMS
    HMMS = []
    for fam in [f for v in FAMILIES.values() for f in v]:
        for d in ("data/raw/hmm_markers", "data/raw/hmm"):
            p = ROOT / d / f"{fam}.hmm"
            if p.exists():
                with pyhmmer.plan7.HMMFile(p) as hf:
                    HMMS.append(next(iter(hf)))
                break


def one(a):
    acc, d = a
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(Path(d) / f"{acc}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by = {_s(s.name): s for s in seqs}
    hits = {}
    for top in pyhmmer.hmmsearch(HMMS, seqs, bit_cutoffs="gathering", cpus=1):
        fam = _s(top.query.accession).split(".")[0]
        for h in top:
            if h.included:
                hits.setdefault(fam, set()).add(_s(h.name))
    row = {"ncbi_assembly": acc}
    for fam in [f for v in FAMILIES.values() for f in v]:
        names = hits.get(fam, set())
        if not names:
            continue
        n = D_ = E_ = K_ = R_ = 0
        for nm in names:
            s = by[nm].textize().sequence.rstrip("*"); n += len(s); D_ += s.count("D"); E_ += s.count("E"); K_ += s.count("K"); R_ += s.count("R")
        row[f"{fam}_acidic_excess"] = (D_ + E_ - K_ - R_) / n; row[f"{fam}_DE"] = (D_ + E_) / n; row[f"{fam}_KR"] = (K_ + R_) / n; row[f"{fam}_nhits"] = len(names)
    return row


if __name__ == "__main__":
    phen, d, out = sys.argv[1:4]; procs = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    accs = [a for a in dict.fromkeys(pd.read_csv(phen).ncbi_assembly.dropna()) if (Path(d) / f"{a}.faa").exists()]
    print(len(accs), "genomes", flush=True)
    with Pool(procs, initializer=init) as p:
        rows = list(p.imap_unordered(one, [(a, d) for a in accs], chunksize=8))
    # DETERMINISM, correction 17, extended. imap_unordered yields a genome as whatever worker finishes it, so
    # row order varied per run and the file was not reproducible from its inputs. Nothing reads these files
    # positionally, so no figure changes; sorting makes the file a function of the inputs too.
    pd.DataFrame(rows).sort_values("ncbi_assembly").reset_index(drop=True).to_csv(out, index=False); print("done")
