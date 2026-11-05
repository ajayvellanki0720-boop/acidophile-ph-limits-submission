"""Usage: python3 scripts/03_system_counts.py annotations.tsv OUT.csv
annotations.tsv columns: genome_id, protein_id, gene, product (genome_id = NCBI assembly)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from acidophile.systems import count_systems, pathway_scores

c = count_systems(pd.read_csv(sys.argv[1], sep="\t"))
c.merge(pathway_scores(c), on="genome_id").to_csv(sys.argv[2], index=False)
print(f"{len(c)} genomes; PROVISIONAL keyword detection, validate with HMMs before claims")
