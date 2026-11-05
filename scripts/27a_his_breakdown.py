"""Per-proteome histidine breakdown (clustered vs isolated, membrane vs soluble) -> DATA_DIR/his_breakdown.csv. Parallel, resumable.
Usage: python3 scripts/27a_his_breakdown.py phenotypes.csv PROTEOME_DIR OUT.csv [processes]"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.proteome import his_breakdown, parse_fasta


def one(a):
    acc, d = a
    return {"ncbi_assembly": acc, **his_breakdown(parse_fasta(Path(d) / f"{acc}.faa"))}


if __name__ == "__main__":
    ph, d, out = sys.argv[1:4]; procs = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    accs = [a for a in dict.fromkeys(pd.read_csv(ph).ncbi_assembly.dropna()) if (Path(d) / f"{a}.faa").exists()]
    old = pd.read_csv(out) if Path(out).exists() else pd.DataFrame(); have = set(old.ncbi_assembly) if len(old) else set()
    todo = [a for a in accs if a not in have]; print(len(accs), "genomes |", len(todo), "to compute", flush=True)
    with Pool(procs) as p:
        rows = list(p.imap_unordered(one, [(a, d) for a in todo], chunksize=4))
    # DETERMINISM, correction 17, extended. imap_unordered yields a genome as whatever worker finishes it, so
    # row order varied per run and the file was not reproducible from its inputs. Nothing reads these files
    # positionally, so no figure changes; sorting makes the file a function of the inputs too.
    pd.concat([old, pd.DataFrame(rows)], ignore_index=True).sort_values("ncbi_assembly").reset_index(drop=True).to_csv(out, index=False); print("done")
