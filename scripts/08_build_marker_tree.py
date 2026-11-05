"""Build per-domain marker trees (+ column-bootstrap trees) for every genome with a downloaded proteome.
Usage: python3 scripts/08_build_marker_tree.py [n_bootstrap] [phenotypes.csv] [outdir]
  defaults: 100, data/curated/phenotypes.csv, results/phylo  -> {Bacteria,Archaea}.nwk, *_bootstrap.nwk, marker_qc.csv"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
import pyhmmer
from acidophile import phylo as ph

B = int(sys.argv[1]) if len(sys.argv) > 1 else 100
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else ROOT / "results/phylo"; OUT.mkdir(parents=True, exist_ok=True)
pheno = pd.read_csv(sys.argv[2] if len(sys.argv) > 2 else ROOT / "data/curated/phenotypes.csv")
have = {p.stem for p in (ROOT / "data/raw/proteomes").glob("*.faa")}
pheno = pheno[pheno.ncbi_assembly.isin(have)].drop_duplicates("ncbi_assembly")
hmms = []
hmm_by = {}
for fam in ph.MARKERS:
    with pyhmmer.plan7.HMMFile(ROOT / f"data/raw/hmm_markers/{fam}.hmm") as hf:
        h = next(iter(hf)); hmms.append(h); hmm_by[fam] = h
hits = {}
for acc in pheno.ncbi_assembly:
    hits[acc] = ph.best_marker_hits(ROOT / f"data/raw/proteomes/{acc}.faa", hmms)
qc = []
for dom in ("Bacteria", "Archaea"):
    genomes = list(pheno[pheno.domain == dom].ncbi_assembly)
    if not genomes:
        continue
    per_family = {}
    for fam in ph.MARKERS:
        seqs = {g: hits[g][fam][1] for g in genomes if fam in hits[g]}
        if len(seqs) >= max(4, 0.5 * len(genomes)):
            per_family[fam] = ph.align_family(hmm_by[fam], seqs)
    M, kept = ph.concatenate(per_family, genomes)
    dropped = sorted(set(genomes) - set(kept))
    D = ph.distance_matrix_fast(M)
    tree = ph.nj_tree_fast(D, kept)
    (OUT / f"{dom}.nwk").write_text(ph.to_newick(tree))
    boots = ph.bootstrap_trees(M, kept, B=B, seed=1)
    (OUT / f"{dom}_bootstrap.nwk").write_text("\n".join(ph.to_newick(t) for t in boots))
    qc.append({"domain": dom, "genomes": len(genomes), "kept": len(kept), "families_used": len(per_family), "alignment_columns": M.shape[1],
               "median_marker_families_per_genome": int(pd.Series([len(hits[g]) for g in kept]).median()), "dropped": ";".join(dropped)})
    print(qc[-1], flush=True)
pd.DataFrame(qc).to_csv(OUT / "marker_qc.csv", index=False)
