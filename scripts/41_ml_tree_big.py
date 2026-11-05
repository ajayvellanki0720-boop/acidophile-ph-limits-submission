"""Maximum-likelihood tree for the full bacterial set (2,768 genomes) -> results/phylo_broad2_ml/.
The biggest standing methodological gap: every result on this dataset rests on neighbour-joining trees of one concatenated ribosomal-protein alignment. The 25 column-bootstrap
trees in 39b show the estimate is stable across replicates of that alignment, but they cannot detect error shared by every neighbour-joining tree of it. This builds the same
29-family concatenation and runs IQ-TREE on it.
Differences from scripts/14_ml_tree.py, which did 235 genomes: the substitution model is FIXED to LG+R4 (what ModelFinder selected by BIC on this same marker set for bacteria)
rather than re-searched, and no bootstrap or SH-aLRT replicates are run. At 2,768 taxa ModelFinder plus 1,000 ultrafast bootstraps would take weeks; a single tree search is
days. The result therefore carries no branch support and is only for re-running the phylogeny-aware tests on a non-distance-based topology.
Stage 1 writes the alignment and prints IQ-TREE's own memory requirement, then stops unless --run is passed, because the run is long enough to be worth checking first.

MEASURED 2026-10-26 on the alignment this builds (2,768 sequences, 3,144 columns, 3,117 distinct patterns): IQ-TREE reports "5367 MB RAM (5 GB) is required". On a machine with
8 GB of RAM that is most of it, for a search expected to take days, so this was NOT launched: it would leave the machine swapping and unusable for other work throughout. Run it
where there is more memory, or accept the cost deliberately. The alignment is written and committed-adjacent so the run can start immediately elsewhere.
Usage: python3 scripts/41_ml_tree_big.py [--run] [threads]
Then:  PHYLO_DIR=results/phylo_broad2_ml python3 scripts/36d_exposure_robustness.py   (and 40c, 39a)"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd
import pyhmmer
from Bio import Phylo
from acidophile import phylo as ph

RUN = "--run" in sys.argv
THREADS = next((a for a in sys.argv[1:] if a.isdigit()), "6")
OUT = ROOT / "results/phylo_broad2_ml"; OUT.mkdir(parents=True, exist_ok=True)
IQ = str(ROOT / "tools/iqtree3")
AFA = OUT / "Bacteria.afa"
PHENO = ROOT / "results/broad2/phenotypes_all.csv"

if not AFA.exists():
    pheno = pd.read_csv(PHENO)
    have = {p.stem for p in (ROOT / "data/raw/proteomes").glob("*.faa")}
    pheno = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(have)].drop_duplicates("ncbi_assembly")
    genomes = list(pheno.ncbi_assembly)
    print(f"{len(genomes)} bacterial genomes; finding marker hits", flush=True)
    hmms, hmm_by = [], {}
    for fam in ph.MARKERS:
        with pyhmmer.plan7.HMMFile(ROOT / f"data/raw/hmm_markers/{fam}.hmm") as hf:
            h = next(iter(hf)); hmms.append(h); hmm_by[fam] = h
    hits = {}
    for i, acc in enumerate(genomes, 1):
        hits[acc] = ph.best_marker_hits(ROOT / f"data/raw/proteomes/{acc}.faa", hmms)
        if i % 250 == 0:
            print(f"  {i}/{len(genomes)}", flush=True)
    per_family = {}
    for fam in ph.MARKERS:
        seqs = {g: hits[g][fam][1] for g in genomes if fam in hits[g]}
        if len(seqs) >= max(4, 0.5 * len(genomes)):
            per_family[fam] = ph.align_family(hmm_by[fam], seqs)
    M, kept = ph.concatenate(per_family, genomes)
    AFA.write_text("".join(f">{n}\n{''.join(row)}\n" for n, row in zip(kept, M)))
    print(f"alignment: {len(kept)} genomes x {M.shape[1]} columns -> {AFA}", flush=True)
else:
    print(f"alignment already present: {AFA}", flush=True)

prefix = OUT / "iq_Bacteria"
cmd = [IQ, "-s", str(AFA), "--prefix", str(prefix), "-m", "LG+R4", "-T", THREADS, "--seed", "20261026", "-redo"]
if not RUN:
    print("\nchecking IQ-TREE's memory requirement (no tree search):", flush=True)
    r = subprocess.run(cmd + ["-n", "0"], capture_output=True, text=True, timeout=3600)
    for line in (r.stdout + r.stderr).splitlines():
        if any(k in line for k in ("Site patterns", "RAM", "memory", "ERROR", "Alignment has")):
            print("  " + line.strip(), flush=True)
    print("\nnot running the tree search; pass --run to start it (expect days).")
    sys.exit(0)

print(f"\nstarting the ML tree search with {THREADS} threads (LG+R4, no bootstrap) - expect days", flush=True)
subprocess.run(cmd, check=True)
t = Phylo.read(str(prefix.with_suffix(".treefile")), "newick"); t.root_at_midpoint()
(OUT / "Bacteria.nwk").write_text(ph.to_newick(t))
print("ALL_DONE: wrote", OUT / "Bacteria.nwk", flush=True)
