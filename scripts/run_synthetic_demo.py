"""End-to-end smoke run on SYNTHETIC data. Results mean nothing biologically; they only show the pipeline
runs and that blocked validation behaves. Real data goes through the same functions."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from acidophile.synthetic import make_synthetic
from acidophile.models import compare_models
from acidophile.lineage import matched_pairs, contrast_features

df = make_synthetic()
res, _ = compare_models(df, scheme="family_blocked")
print(res.round(3).to_string(index=False))
pairs = matched_pairs(df, rank="genus")
print(f"\nmatched pairs: {len(pairs)}")
print(contrast_features(df, pairs, ["mean_pI", "acidic_excess", "rate_kdp", "pw_membrane_barrier"], min_pairs=4).round(3).to_string(index=False))
