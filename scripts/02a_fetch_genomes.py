"""Download genome + protein FASTA for each accession in a phenotype CSV (NCBI Datasets, public API).
Files land in data/raw/genomes/<acc>.fna and data/raw/proteomes/<acc>.faa after md5 verification against the
archive's md5sum.txt. Also writes data/curated/annotations.tsv (protein_id, product) parsed from FASTA headers.
Usage: python3 scripts/02a_fetch_genomes.py data/curated/phenotypes_primary.csv [annotations_out.tsv]"""
import hashlib
import io
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GEN, PROT = ROOT / "data/raw/genomes", ROOT / "data/raw/proteomes"
URL = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download"
       "?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&filename=x.zip")


def fetch(acc):
    z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(URL.format(acc), timeout=300).read()))
    sums = {l.split()[1].lstrip("*"): l.split()[0] for l in z.read("md5sum.txt").decode().splitlines() if l.strip()}
    out = {}
    for name in z.namelist():
        if name.endswith("_genomic.fna") or name.endswith("protein.faa"):
            data = z.read(name)
            if hashlib.md5(data).hexdigest() != sums.get(name):
                raise ValueError(f"md5 mismatch for {name}")
            out["fna" if name.endswith(".fna") else "faa"] = data
    return out


if __name__ == "__main__":
    GEN.mkdir(parents=True, exist_ok=True); PROT.mkdir(parents=True, exist_ok=True)
    accs = pd.read_csv(sys.argv[1])["ncbi_assembly"].dropna().unique()
    ann, failed, used_log = [], [], []
    for a in accs:
        f, p = GEN / f"{a}.fna", PROT / f"{a}.faa"
        if not (f.exists() and p.exists()):
            try:
                d, used = fetch(a), a
                if ("faa" not in d or "fna" not in d) and a.startswith("GCA_"):
                    used = "GCF_" + a[4:]          # RefSeq copy of the same assembly carries the annotation
                    d = fetch(used)
                if "faa" not in d or "fna" not in d:
                    raise ValueError("no protein.faa (assembly not annotated)")
                f.write_bytes(d["fna"]); p.write_bytes(d["faa"])
                used_log.append((a, used))
            except Exception as e:
                failed.append((a, str(e))); print(a, "FAILED", e, flush=True); continue
            time.sleep(0.4)
        for line in p.read_text().splitlines():
            if line.startswith(">"):
                pid, _, rest = line[1:].partition(" ")
                ann.append((a, pid, "", rest.split(" [")[0]))
        print(a, "ok", flush=True)
    pd.DataFrame(used_log, columns=["requested", "fetched"]).to_csv(ROOT / "data/curated/accession_used.csv", mode="a", index=False, header=False)
    pd.DataFrame(ann, columns=["genome_id", "protein_id", "gene", "product"]).to_csv(Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data/curated/annotations.tsv", sep="\t", index=False)
    print(f"{len(accs) - len(failed)} downloaded, {len(failed)} failed: {failed}")
