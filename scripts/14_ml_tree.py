"""Maximum-likelihood marker trees with IQ-TREE 3.1.4 (official macOS arm64 release, SHA-256 verified against the digest GitHub
publishes for the asset; binary at tools/iqtree3, git-ignored).

Same concatenated ribosomal-protein alignments as the neighbour-joining trees (29 Pfam families, gathering-cutoff best hits).
Model selection: ModelFinder (-m MFP) over LG, WAG, JTT and Q.pfam with and without empirical frequencies and with gamma / FreeRate
heterogeneity chosen by BIC; 1000 ultrafast bootstraps (-B) and 1000 SH-aLRT replicates; UFBoot trees are kept with branch lengths
(-wbtl) so downstream phylogenetic tests can be repeated over them. Midpoint rooted for the regressions.

Usage: python3 scripts/14_ml_tree.py [phenotypes.csv] [outdir] [threads]
Then:  PHYLO_DIR=results/phylo_ml python3 scripts/12_full_phylo_screen.py   (and 13_kdp_trk_diagnostics.py, 15_tree_sensitivity.py)"""
import io
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
import pyhmmer
from Bio import Phylo
from acidophile import phylo as ph

PHENO = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/full/phenotypes_all.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results/phylo_ml"
THREADS = sys.argv[3] if len(sys.argv) > 3 else "6"
IQ = str(ROOT / "tools/iqtree3")
OUT.mkdir(parents=True, exist_ok=True)


def to_fasta(M, names):
    return "".join(f">{n}\n{''.join(row)}\n" for n, row in zip(names, M))


def midpoint(path: Path):
    t = Phylo.read(str(path), "newick"); t.root_at_midpoint(); return t


pheno = pd.read_csv(PHENO)
have = {p.stem for p in (ROOT / "data/raw/proteomes").glob("*.faa")}
pheno = pheno[pheno.ncbi_assembly.isin(have)].drop_duplicates("ncbi_assembly")
hmms, hmm_by = [], {}
for fam in ph.MARKERS:
    with pyhmmer.plan7.HMMFile(ROOT / f"data/raw/hmm_markers/{fam}.hmm") as hf:
        h = next(iter(hf)); hmms.append(h); hmm_by[fam] = h
hits = {a: ph.best_marker_hits(ROOT / f"data/raw/proteomes/{a}.faa", hmms) for a in pheno.ncbi_assembly}
for dom in ("Bacteria", "Archaea"):
    genomes = list(pheno[pheno.domain == dom].ncbi_assembly)
    per_family = {}
    for fam in ph.MARKERS:
        seqs = {g: hits[g][fam][1] for g in genomes if fam in hits[g]}
        if len(seqs) >= max(4, 0.5 * len(genomes)):
            per_family[fam] = ph.align_family(hmm_by[fam], seqs)
    M, kept = ph.concatenate(per_family, genomes)
    afa = OUT / f"{dom}.afa"
    afa.write_text(to_fasta(M, kept))
    print(f"{dom}: {len(kept)} genomes x {M.shape[1]} columns -> running IQ-TREE", flush=True)
    prefix = OUT / f"iq_{dom}"
    subprocess.run([IQ, "-s", str(afa), "--prefix", str(prefix), "-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam", "-mrate", "E,G,R4",
                    "-B", "1000", "-alrt", "1000", "-wbtl", "-T", THREADS, "--seed", "20261024", "-redo", "--quiet"], check=True)
    main = midpoint(prefix.with_suffix(".treefile"))
    (OUT / f"{dom}.nwk").write_text(ph.to_newick(main))
    boots = list(Phylo.parse(str(prefix.with_suffix(".ufboot")), "newick"))
    step = max(1, len(boots) // 100)
    reps = []
    for t in boots[::step][:100]:
        t.root_at_midpoint(); reps.append(ph.to_newick(t))
    (OUT / f"{dom}_bootstrap.nwk").write_text("\n".join(reps))
    best = [l for l in prefix.with_suffix(".iqtree").read_text().splitlines() if l.startswith(("Best-fit model", "Log-likelihood of the tree"))]
    print(f"{dom}: done | " + " | ".join(best) + f" | {len(reps)} UFBoot trees kept", flush=True)
print("ALL_DONE", flush=True)
