"""How sensitive is the Bacteria Kdp / Trk-Ktr PGLS result to how the tree was built? Uses only in-Python tree methods: -> DATA_DIR/tree_sensitivity.csv
NJ with other gamma shapes, UPGMA, trees from independent halves of the marker set (large-subunit vs small-subunit ribosomal
proteins), and the 100 column-bootstrap NJ trees. Usage: python3 scripts/15_tree_sensitivity.py"""
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import pyhmmer
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D = ROOT / "results/full"
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tips_all = [t.name for t in Phylo.read(ROOT / "results/phylo_full/Bacteria.nwk", "newick").get_terminals()]
df = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(tips_all)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
cov = ["ph_opt", "src_thermobase", "src_acidb", "src_bacdive"]
y = df.ph_min.to_numpy(float)
z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)
Xbase = np.column_stack([np.ones(len(df))] + [df[c].to_numpy(float) for c in cov])
ZK, ZT = z(df.hm_rate_kdp), z(df.hm_rate_trk_ktr)


def fit(tree):
    C = ph.vcv_fast(tree, list(df.ncbi_assembly))
    j = ph.pgls(y, np.column_stack([Xbase, ZK, ZT]), C)
    k = Xbase.shape[1]
    return {"lambda": j["lambda"], "kdp_b": j["beta"][k], "kdp_p": j["p"][k], "trk_b": j["beta"][k + 1], "trk_p": j["p"][k + 1]}


rows = []


def add(label, tree):
    r = fit(tree); r["tree"] = label; rows.append(r)


# markers -> alignment (same procedure as the main tree)
pheno = df[["ncbi_assembly"]]
hmms, hmm_by = [], {}
for fam in ph.MARKERS:
    with pyhmmer.plan7.HMMFile(ROOT / f"data/raw/hmm_markers/{fam}.hmm") as hf:
        h = next(iter(hf)); hmms.append(h); hmm_by[fam] = h
hits = {a: ph.best_marker_hits(ROOT / f"data/raw/proteomes/{a}.faa", hmms) for a in df.ncbi_assembly}
genomes = list(df.ncbi_assembly)
per_family = {}
for fam in ph.MARKERS:
    seqs = {g: hits[g][fam][1] for g in genomes if fam in hits[g]}
    if len(seqs) >= 0.5 * len(genomes):
        per_family[fam] = ph.align_family(hmm_by[fam], seqs)
LSU = "PF00181 PF03947 PF00297 PF00573 PF00673 PF00347 PF00238 PF00252 PF00861 PF00237 PF00276 PF00828 PF00572 PF00298 PF00687 PF00466".split()
SSU = [f for f in ph.MARKERS if f not in LSU]


def matrix(fams):
    M, kept = ph.concatenate({f: per_family[f] for f in fams if f in per_family}, genomes)
    assert kept == genomes, "genome dropped; ordering would change"
    return M


Mall = matrix(ph.MARKERS)
add("main NJ tree (gamma shape 1)", Phylo.read(ROOT / "results/phylo_full/Bacteria.nwk", "newick"))
for a in (0.5, 2.0):
    add(f"NJ, gamma shape {a}", ph.nj_tree(ph.distance_matrix(Mall, alpha=a), genomes))
D1 = ph.distance_matrix(Mall)
lower = [list(D1[i, :i + 1]) for i in range(len(genomes))]
add("UPGMA (gamma shape 1)", DistanceTreeConstructor().upgma(DistanceMatrix(genomes, lower)))
add("NJ, large-subunit ribosomal proteins only", ph.nj_tree(ph.distance_matrix(matrix(LSU)), genomes))
add("NJ, small-subunit ribosomal proteins only", ph.nj_tree(ph.distance_matrix(matrix(SSU)), genomes))
boots = list(Phylo.parse(ROOT / "results/phylo_full/Bacteria_bootstrap.nwk", "newick"))
B = pd.DataFrame([fit(b) for b in boots])
R = pd.DataFrame(rows)[["tree", "lambda", "kdp_b", "kdp_p", "trk_b", "trk_p"]]
R.to_csv(D / "tree_sensitivity.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))
print(f"\n{len(B)} bootstrap NJ trees: Kdp beta {B.kdp_b.quantile(.025):.3f} to {B.kdp_b.quantile(.975):.3f}, p median {B.kdp_p.median():.1e}, max {B.kdp_p.max():.1e}, "
      f"share p<0.05 {(B.kdp_p < .05).mean():.2f} | Trk/Ktr beta {B.trk_b.quantile(.025):.3f} to {B.trk_b.quantile(.975):.3f}, p median {B.trk_p.median():.1e}, "
      f"max {B.trk_p.max():.1e}, share p<0.05 {(B.trk_p < .05).mean():.2f} | lambda {B['lambda'].min():.2f}-{B['lambda'].max():.2f}")
B.to_csv(D / "tree_sensitivity_bootstrap.csv", index=False)
