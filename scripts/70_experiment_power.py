"""Power and decision thresholds for the charge-reversal experiment, fixed before any data exist -> DATA_DIR/experiment_power.csv
docs/CHARGE_REVERSAL_DESIGN.md specifies what to build. It does not fix the analysis, and this project has recorded twelve corrections of its own, most
of which were post-hoc flexibility of exactly the kind a pre-registration exists to prevent: a scale chosen after seeing the answer, a subset selected on
the outcome, a ratio computed across two samples. An experiment designed by the same process, analysed with the same freedom, would inherit the problem.
So the numbers that decide the experiment are computed here, before the bench work, and docs/PREREGISTRATION.md commits to them.
Three things this makes explicit that the design document leaves open.
  1. WHICH PREDICTED EFFECT THE STUDY IS POWERED FOR. The design quotes -0.16 pH units for the primary variant, scaled through the 7-reference panel,
     and notes that the same calculation on the 32-reference panel gives -0.05. That threefold gap is unresolved and the comparative data do not choose
     between the two (docs/NUMBERS.md: the two panels are different measurements and their per-unit figures are not comparable). Powering for the larger
     and observing nothing produces an uninterpretable null, which is the worst outcome available, so the required sample size is computed for BOTH and
     the pre-registration commits to the conservative one.
  2. MULTIPLICITY ACROSS ARMS. The design declares one primary variant, a nine-point dose series, six controls and a reverse arm. That is roughly
     seventeen comparisons at alpha 0.05 each. scripts/68 showed that correction is not what distinguishes a real result from a false one in this
     project's history, but it is still the difference between one pre-specified test and seventeen chances, so the primary endpoint is designated as a
     single test and everything else is explicitly secondary, with the correction stated rather than applied silently later.
  3. WHAT THE POOLED ASSAY CANNOT YET BE POWERED FOR. The primary readout is relative fitness in a barcoded pool, not pH_min. Converting a fitness
     difference into a shift in the growth limit needs a slope - fitness per pH unit of limit - that no one has measured, and inventing one would be the
     same error as the effect sizes this project has already had to withdraw. So the pooled assay's power is left as a function of that unknown and the
     pilot that must measure it is specified; only the confirmation assay, whose readout is pH_min from a fitted growth curve, is powered in pH units.
Usage: python3 scripts/70_experiment_power.py   (env DATA_DIR=results/broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
# the two candidate predicted effects, taken from the design document's own scaling of the two panels
EFFECTS = {"32-reference scaling (conservative)": 0.05, "7-reference scaling (design's primary figure)": 0.16}
SD_ASSAY = 0.10          # SD of pH_min from a fitted growth curve, the design document's stated figure
POWER = 0.80
# the declared comparison structure, which is what sets the multiplicity
ARMS = {"primary combined variant": 1, "dose series": 9, "controls": 6, "reverse arm": 1}
rows = []


def n_per_group(delta, sd, alpha, power=POWER):
    """Two-sample two-sided t-test sample size, solved by search on the exact non-central t power rather than the
    normal approximation, which understates n at the sample sizes in play here."""
    for n in range(3, 100000):
        df = 2 * n - 2
        ncp = delta / (sd * np.sqrt(2.0 / n))
        crit = stats.t.ppf(1 - alpha / 2, df)
        pw = stats.nct.sf(crit, df, ncp) + stats.nct.cdf(-crit, df, ncp)
        if pw >= power:
            return n, float(pw)
    return np.nan, np.nan


def mde(n, sd, alpha, power=POWER):
    """Smallest effect detectable at this n, by bisection on the same power function."""
    lo, hi = 1e-4, 5.0
    for _ in range(80):
        mid = (lo + hi) / 2
        df = 2 * n - 2
        ncp = mid / (sd * np.sqrt(2.0 / n))
        crit = stats.t.ppf(1 - alpha / 2, df)
        pw = stats.nct.sf(crit, df, ncp) + stats.nct.cdf(-crit, df, ncp)
        if pw < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


total_arms = sum(ARMS.values())
print(f"the declared comparison structure: " + ", ".join(f"{k} x{v}" for k, v in ARMS.items()) + f" = {total_arms} comparisons")
ALPHAS = {"0.05, one pre-specified primary test": 0.05,
          f"0.05/{total_arms} Bonferroni across every declared arm": 0.05 / total_arms}

print(f"\n1. required replicates per group, confirmation assay (readout pH_min from a fitted curve, SD {SD_ASSAY})")
print(f"   {'predicted effect':46s} {'alpha basis':44s} {'n/group':>8s}")
for ename, delta in EFFECTS.items():
    for aname, alpha in ALPHAS.items():
        n, pw = n_per_group(delta, SD_ASSAY, alpha)
        print(f"   {ename:46s} {aname:44s} {n:>8d}")
        # stored SIGNED, because the prediction is a decrease in pH_min and an unsigned 0.05 made the
        # pre-registration's "-0.05" untraceable to this file - caught by scripts/69 on that document's first run
        rows.append({"part": "required n", "effect": ename, "delta_ph": -delta, "alpha_basis": aname, "alpha": alpha,
                     "sd_assay": SD_ASSAY, "n_per_group": int(n), "achieved_power": pw})

print(f"\n2. what the design's stated floor of n>=20 actually buys")
for aname, alpha in ALPHAS.items():
    m = mde(20, SD_ASSAY, alpha)
    print(f"   at n=20 per group, alpha {aname:44s} smallest detectable shift {m:.3f} pH units")
    rows.append({"part": "mde at design n", "n_per_group": 20, "alpha_basis": aname, "alpha": alpha,
                 "sd_assay": SD_ASSAY, "mde_ph": float(m)})
m20 = mde(20, SD_ASSAY, 0.05)
small = min(EFFECTS.values())
print(f"   The design's floor was set against the -0.16 prediction. Against the conservative -{small:.2f} it is "
      f"{'ADEQUATE' if m20 <= small else 'NOT ADEQUATE'}:")
print(f"   n=20 detects {m20:.3f} and the conservative prediction is {small:.2f}, so the study would be underpowered by a factor of "
      f"{m20 / small:.1f} on that scaling.")
n_cons, _ = n_per_group(small, SD_ASSAY, 0.05)
print(f"   Powering for the conservative prediction instead needs n={n_cons} per group at alpha 0.05, and "
      f"n={n_per_group(small, SD_ASSAY, ALPHAS[list(ALPHAS)[1]])[0]} with the arm correction.")
rows.append({"part": "verdict", "note": "design floor n=20 against the conservative prediction",
             "mde_at_20": float(m20), "conservative_delta": -small, "shortfall_factor": float(m20 / small),
             "n_required_conservative": int(n_cons)})

print(f"\n3. sensitivity of that verdict to the assay SD, which is itself a pilot measurement")
print(f"   {'assay SD':>10s} {'n for -0.05':>12s} {'n for -0.16':>12s} {'MDE at n=20':>13s}")
for sd in (0.05, 0.075, 0.10, 0.15, 0.20):
    n_small, _ = n_per_group(0.05, sd, 0.05)
    n_big, _ = n_per_group(0.16, sd, 0.05)
    print(f"   {sd:>10.3f} {n_small:>12d} {n_big:>12d} {mde(20, sd, 0.05):>13.3f}")
    rows.append({"part": "sd sensitivity", "sd_assay": sd, "delta_small": -0.05, "delta_large": -0.16,
                 "n_for_0.05": int(n_small), "n_for_0.16": int(n_big),
                 "mde_at_20": float(mde(20, sd, 0.05))})
print("   The design states SD 0.1 without a source. It is the single parameter that decides whether this experiment is")
print("   feasible at a reasonable scale, so the pilot must measure it before the panel is built, not after.")

# the PER-REVERSAL prediction, which the design document quotes as -0.011 and derives in prose. A derived figure
# quoted anywhere has to be a stored row or it is unguarded (correction 12), so it is computed here from the
# committed combined-variant prediction and the 14 reversed positions.
N_POSITIONS = 14
print("\n3b. the per-reversal prediction implied by each scaling, over the 14 reversed positions")
for ename, delta in EFFECTS.items():
    per = -delta / N_POSITIONS
    print(f"   {ename:46s} {per:+.4f} pH units per reversal")
    rows.append({"part": "per reversal", "effect": ename, "delta_ph": -delta, "n_positions": N_POSITIONS,
                 "per_reversal_ph": float(per)})
print("   The design document quotes -0.011, which is the 7-reference figure. On the conservative scaling it is -0.0036,")
print("   an order of magnitude below any assay resolution yet demonstrated, which is why the dose series is secondary.")

print("\n4. the pooled assay cannot be powered yet, and this says exactly what is missing")
print("   Its readout is log relative fitness, not pH_min. Converting between them needs a slope - fitness per pH unit of")
print("   growth limit - that has not been measured for this chassis or any other. Rather than invent one:")
print("   REQUIRED PILOT. Compete the reconstructed wild type against a small ladder of strains whose pH_min differs by")
print("   independently measured amounts, across the same pH series, and regress log fitness on known pH_min difference.")
print("   That slope, with its standard error, is what makes the pooled assay's resolution a number instead of a hope.")
print("   Until it exists, the design's claim to resolve 0.02 pH units is an assumption, and the pre-registration records it")
print("   as the gate the experiment must pass before the dose series is interpreted at all.")
rows.append({"part": "pooled assay", "note": "not powerable without a measured fitness-per-pH-unit slope",
             "required_pilot": "regress log relative fitness on independently measured pH_min differences across a spike-in ladder"})

R = pd.DataFrame(rows)
R.to_csv(D / "experiment_power.csv", index=False)
print("\nwrote", D / "experiment_power.csv")
print("\nReading. The headline number is the shortfall in part 2: the design's replicate floor was set against the larger of two")
print("predicted effects that the comparative data do not choose between, and against the smaller one it is underpowered several")
print("fold. That is a decision to take before building rather than a result to explain afterwards, which is the whole point of")
print("computing it now.")
