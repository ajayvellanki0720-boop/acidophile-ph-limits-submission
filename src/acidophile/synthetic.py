"""Synthetic phenotype+genome table for tests and smoke runs. NOT data."""
import numpy as np
import pandas as pd


def make_synthetic(n_families=14, per_family=6, seed=0, signal=True):
    rng = np.random.default_rng(seed)
    rows = []
    for f in range(n_families):
        fam_shift = rng.normal(0, 0.4)
        for k in range(per_family):
            opt = np.clip(rng.normal(3.0, 1.2), 0.8, 6.5)
            kdp = rng.poisson(1 + (2 if f % 3 == 0 else 0))
            mean_pi = rng.normal(6.0, 0.6)
            ph_min = opt - 1.0 - (0.35 * kdp if signal else 0) + fam_shift + rng.normal(0, 0.2)
            rows.append({
                "strain_id": f"S{f:02d}_{k}", "family": f"Fam{f}", "genus": f"Gen{f}_{k % 3}",
                "class": f"Cls{f % 5}", "domain": "Bacteria" if f % 4 else "Archaea",
                "acid_identity": "H2SO4" if f % 2 else "HCl", "temp_opt_c": rng.normal(40, 15),
                "ph_opt": opt, "ph_min": ph_min, "ph_max": opt + 2, "ph_min_bound": "measured",
                "confidence": "high", "ncbi_assembly": f"GCA_{f:03d}{k:03d}.1",
                "mean_pI": mean_pi, "acidic_excess": rng.normal(0, 0.02), "gc": rng.uniform(0.3, 0.7),
                "rate_kdp": float(kdp), "rate_gad": float(rng.poisson(1)), "hm_rate_kdp": float(kdp), "hm_rate_gad": float(rng.poisson(1)),
                "pw_membrane_barrier": float(rng.integers(0, 2)), "hm_pw_membrane_barrier": float(rng.integers(0, 2)),
            })
    return pd.DataFrame(rows)
