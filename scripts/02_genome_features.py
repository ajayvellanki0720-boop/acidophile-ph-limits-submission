"""Proteome + genome features for each assembly (parallel; reuses rows already in OUT.csv).
Usage: python3 scripts/02_genome_features.py phenotypes.csv PROTEOME_DIR GENOME_DIR OUT.csv [processes]
Expects PROTEOME_DIR/<assembly>.faa and GENOME_DIR/<assembly>.fna (fetch with scripts/02a_fetch_genomes.py)."""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.proteome import genome_features, parse_fasta, proteome_features


def one(args):
    acc, prot_dir, gen_dir = args
    faa, fna = Path(prot_dir) / f"{acc}.faa", Path(gen_dir) / f"{acc}.fna"
    if not (faa.exists() and fna.exists()):
        return acc, None
    return acc, {"ncbi_assembly": acc, **proteome_features(parse_fasta(faa)), **genome_features(fna)}


if __name__ == "__main__":
    ph, prot_dir, gen_dir, out = sys.argv[1:5]
    procs = int(sys.argv[5]) if len(sys.argv) > 5 else 6
    accs = list(dict.fromkeys(pd.read_csv(ph)["ncbi_assembly"].dropna()))
    old = pd.read_csv(out) if Path(out).exists() else pd.DataFrame()
    have = set(old.ncbi_assembly) if len(old) else set()
    todo = [a for a in accs if a not in have]
    print(f"{len(accs)} assemblies requested | {len(have & set(accs))} already computed | {len(todo)} to compute with {procs} processes", flush=True)
    rows, missing = [], []
    with Pool(procs) as pool:
        for i, (acc, r) in enumerate(pool.imap_unordered(one, [(a, prot_dir, gen_dir) for a in todo], chunksize=4), 1):
            (rows.append(r) if r else missing.append(acc))
            if i % 100 == 0:
                print(f"  {i}/{len(todo)}", flush=True)
    res = pd.concat(([old[old.ncbi_assembly.isin(accs)]] if len(old) else []) + [pd.DataFrame(rows)], ignore_index=True)
    # DETERMINISM, correction 17, extended. imap_unordered yields a genome as whatever worker finishes it, so
    # row order varied per run and the file was not reproducible from its inputs. Nothing reads these files
    # positionally, so no figure changes; sorting makes the file a function of the inputs too.
    res = res.sort_values("ncbi_assembly").reset_index(drop=True)
    res.to_csv(out, index=False)
    print(f"{len(res)} assemblies in {out}; {len(missing)} missing files: {missing[:10]}")
