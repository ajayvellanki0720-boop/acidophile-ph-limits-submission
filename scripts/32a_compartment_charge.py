"""Per-proteome charge composition by protein class (membrane / secreted-like / cytoplasmic / conserved core / non-core soluble) -> DATA_DIR/compartment_charge.csv.
Conserved core = the best hit of each of the 29 ribosomal-protein Pfam markers. Parallel, resumable.
Usage: python3 scripts/32a_compartment_charge.py phenotypes.csv PROTEOME_DIR OUT.csv [processes]"""
import sys
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
import pyhmmer
from acidophile import phylo as ph
from acidophile.proteome import compartment_composition, parse_fasta

HMMS = None


def init():
    global HMMS
    HMMS = []
    for fam in ph.MARKERS:
        with pyhmmer.plan7.HMMFile(ROOT / f"data/raw/hmm_markers/{fam}.hmm") as hf:
            HMMS.append(next(iter(hf)))


def one(a):
    acc, d = a
    path = Path(d) / f"{acc}.faa"
    names = {n for n, _ in ph.best_marker_hits(path, HMMS).values()}
    return {"ncbi_assembly": acc, **compartment_composition(parse_fasta(path), names)}


if __name__ == "__main__":
    phen, d, out = sys.argv[1:4]; procs = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    accs = [a for a in dict.fromkeys(pd.read_csv(phen).ncbi_assembly.dropna()) if (Path(d) / f"{a}.faa").exists()]
    old = pd.read_csv(out) if Path(out).exists() else pd.DataFrame(); have = set(old.ncbi_assembly) if len(old) else set()
    todo = [a for a in accs if a not in have]; print(len(accs), "genomes |", len(todo), "to compute", flush=True)
    with Pool(procs, initializer=init) as p:
        rows = []
        for i, r in enumerate(p.imap_unordered(one, [(a, d) for a in todo], chunksize=4), 1):
            rows.append(r)
            if i % 250 == 0:
                print(f"  {i}/{len(todo)}", flush=True)
    # DETERMINISM, correction 17, extended. imap_unordered yields a genome as whatever worker finishes it, so
    # row order varied per run and the file was not reproducible from its inputs. Nothing reads these files
    # positionally, so no figure changes; sorting makes the file a function of the inputs too.
    pd.concat([old, pd.DataFrame(rows)], ignore_index=True).sort_values("ncbi_assembly").reset_index(drop=True).to_csv(out, index=False); print("done")
