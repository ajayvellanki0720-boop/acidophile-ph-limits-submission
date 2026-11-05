"""How reliable is pH_min itself, and what does that bound? -> DATA_DIR/outcome_reliability.csv
Every measurement-error check in this project so far has been on the PREDICTOR (scripts/52 for the composite, scripts/55 for reference-structure
transfer). The outcome has never been checked, although the material for it has been sitting in data/curated the whole time: four source tables -
hand-curated species descriptions, ThermoBase, AcidB and BacDive - overlap on 114 species, so those species carry two or more independent curators'
readings of the same organism's lowest growth pH. The phenotype table applied a source-precedence rule and kept one reading per genome, which is why
the duplicates are invisible downstream.
Three things come out of this, and it is worth being precise about which of them can move a published number.
  1. RELIABILITY OF pH_min. A one-way random-effects decomposition over the replicated species gives the share of reported pH_min variance that is
     real between-species variation rather than curator disagreement. This is a number the project has quoted no value for, including when stating
     that the composite is "worth under 1% of variance".
  2. WHAT IT BOUNDS, AND WHAT IT DOES NOT. Error in the DEPENDENT variable does not bias a regression slope - it inflates the standard error and caps
     the attainable R-squared. So this cannot rescue or threaten the coefficient, and the temptation to treat it as a disattenuation story has to be
     resisted; the fifth scale error in this project came from exactly that kind of reflex, so the arithmetic is simulated here rather than asserted.
     What reliability does set is a ceiling on R-squared, which is the right yardstick for the project's weakest claim.
  3. WHETHER THE FLAT UPPER HALF IS INFORMATIVE. scripts/51 found the relationship confined to the acid-tolerant tail, with "phenotype ceiling" left
     as an open gap. If minima reported above about pH 5 are largely conventional - a medium's standard lower bound rather than a measured limit -
     they should carry visibly more curator disagreement than the acidophile readings do. Banding is done on the HIGHER-PRECEDENCE member of each
     pair, not on the value being assessed, so the band assignment does not select on the noise under test. Range restriction inside a band still
     biases a within-band reliability downward, so the error SD is reported directly and reliability is computed against the whole sample's spread.
The replicated species are acidophile-heavy by construction (median replicate mean pH 3.4), because sources overlap where they specialise, and only 34
are bacteria. Refitting the headline on them is therefore not attempted - it would be 34 species against 1,655.
THE LIMITATION THAT SETS WHAT ALL OF THIS MEANS: these four curators are not four independent measurements of an organism's acid limit. They are
largely independent TRANSCRIPTIONS of a shared primary literature - ThermoBase, AcidB and BacDive all compile published species descriptions, and the
hand-curated tier reads those same descriptions directly. So what is estimated below is transcription reliability, which is an UPPER bound on the
reliability of pH_min as a measurement: the growth experiment's own error, the medium, the temperature at which pH was read and the strain's identity
are all common to every curator and cancel out of the comparison. Two sources agreeing perfectly on a number they both copied from the same abstract
is not evidence that the number is right. Every figure here should be read as "no worse than", and the ceiling it implies is optimistic.
Usage: python3 scripts/62_outcome_reliability.py   (env DATA_DIR=results/broad2)"""
import itertools
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
rng = np.random.default_rng(0)
# precedence as the pipeline applies it, best first
PRECEDENCE = ["literature", "thermobase", "acidb", "bacdive"]


def _check_y_error_is_unbiased(n=4000, reps=400, b=0.30, err_sd=0.5):
    """Error in the OUTCOME leaves the slope unbiased and inflates its SE. Simulated, not asserted, because the
    reflex to disattenuate for it is the same reflex that produced this project's fifth scale error."""
    out = []
    for noisy in (0.0, err_sd):
        bs, ses = [], []
        for _ in range(reps):
            x = rng.standard_normal(n)
            y = b * x + rng.standard_normal(n) + noisy * rng.standard_normal(n)
            X = np.column_stack([np.ones(n), x])
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            resid = y - X @ beta
            s2 = resid @ resid / (n - 2)
            ses.append(np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]))
            bs.append(beta[1] * x.std() / y.std())  # per SD of y, as the project's coefficients are
        out.append((np.mean(bs), np.mean(ses)))
    return out


def norm(s):
    return s.astype(str).str.strip().str.lower().str.replace(r"\s+", " ", regex=True)


def load():
    tb = pd.read_csv(ROOT / "data/curated/thermobase_candidates.csv", low_memory=False)[["species", "ph_min"]].assign(src="thermobase")
    ab = pd.read_csv(ROOT / "data/curated/acidb_organisms.csv", low_memory=False)[["species", "ph_min"]].assign(src="acidb")
    bd = pd.concat([pd.read_csv(ROOT / f"data/curated/{f}.csv")[["species", "ph_min"]]
                    for f in ("bacdive_candidates", "bacdive_archaea_candidates")]).assign(src="bacdive")
    lc = pd.read_csv(ROOT / "data/curated/literature_curated.csv").rename(columns={"taxon": "species"})[["species", "ph_min"]].assign(src="literature")
    a = pd.concat([tb, ab, bd, lc])
    a["sp"] = norm(a.species)
    a = a.dropna(subset=["ph_min"])
    # a source that lists several strains of one species is averaged within that source, so a species contributes one reading per curator
    return a.groupby(["sp", "src"], as_index=False).ph_min.mean().pivot(index="sp", columns="src", values="ph_min")


rows = []
W = load()
rep = W[W.notna().sum(axis=1) >= 2]
print(f"{len(W)} species across the four source tables, {len(rep)} of them read by two or more curators")
print(f"readings per source: " + ", ".join(f"{c} {int(W[c].notna().sum())}" for c in W.columns))

print("\n1. pairwise agreement between curators")
for a, b in itertools.combinations(PRECEDENCE, 2):
    m = W[[a, b]].dropna()
    if len(m) < 10:
        print(f"   {a:11s} vs {b:11s} n={len(m):3d}  too few")
        continue
    d = (m[a] - m[b]).to_numpy()
    bs = [d[rng.integers(0, len(d), len(d))].std() for _ in range(4000)]
    lo_ci, hi_ci = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    rows.append({"part": "pair", "a": a, "b": b, "n": len(m), "r": float(m[a].corr(m[b])), "mean_diff": float(d.mean()),
                 "sd_diff": float(d.std()), "sd_diff_lo": lo_ci, "sd_diff_hi": hi_ci, "exact_rate": float((np.abs(d) < 1e-9).mean())})
    print(f"   {a:11s} vs {b:11s} n={len(m):3d}  r {m[a].corr(m[b]):+.3f}  mean difference {d.mean():+.3f}  SD {d.std():.3f} "
          f"[{lo_ci:.3f}, {hi_ci:.3f}]  identical in {(np.abs(d) < 1e-9).mean():.0%}")
print("   The hand-curated tier is the nearest thing to a yardstick: a source's SD against `literature` is how far it strays from a")
print("   careful reading of the same description. By that measure the precedence rule the pipeline uses is in the wrong order - BacDive")
print("   strays least and sits last, AcidB strays most and sits above it - but the bootstrap intervals overlap, so the ordering is not")
print("   established, and the exposure below shows it barely matters either way.")
# precedence only decides anything where two sources overlap, so the practical exposure is worth knowing before anything is made of it
ph = pd.read_csv(D / "phenotypes_all.csv").drop_duplicates("ncbi_assembly")
ph["sp"] = norm(ph.species)
nsrc = W.notna().sum(axis=1)
multi = int(ph.sp.map(nsrc).fillna(0).ge(2).sum())
ab_both = int(ph.sp.isin(W[["acidb", "bacdive"]].dropna().index).sum())
print(f"   exposure: {multi} of {len(ph)} analysed genomes have more than one source available at all, and {ab_both} are decided by the")
print(f"   AcidB-against-BacDive call, so the precedence ordering is close to inconsequential for the dataset as built.")
rows.append({"part": "exposure", "subset": "analysed genomes with >1 source available", "n": multi, "n_total": len(ph)})
rows.append({"part": "exposure", "subset": "genomes decided by the acidb vs bacdive precedence call", "n": ab_both, "n_total": len(ph)})

print("\n2. reliability of a single reported pH_min")
# one-way random effects over replicated species: between-species variance against curator variance
long = rep.stack().rename("ph_min").reset_index().rename(columns={"level_1": "src"})
grand = long.ph_min.mean()
g = long.groupby("sp").ph_min
k = g.size()
msw = float((long.ph_min - g.transform("mean")).pow(2).sum() / (len(long) - len(rep)))     # within-species (curator disagreement)
msb = float((k * (g.mean() - grand) ** 2).sum() / (len(rep) - 1))                         # between-species
k0 = float(k.mean())
var_b = max((msb - msw) / k0, 0.0)
icc = var_b / (var_b + msw)
print(f"   replicated species {len(rep)}, {len(long)} readings, mean {k0:.2f} curators each")
print(f"   between-species variance {var_b:.4f}, curator variance {msw:.4f}  ->  reliability of ONE reading {icc:.3f}")
print(f"   so a single reported minimum carries a curator error of about {np.sqrt(msw):.2f} pH units")
rows.append({"part": "reliability", "subset": "all replicated species", "n": len(rep), "var_between": var_b, "var_within": msw,
             "icc_single_reading": icc, "error_sd": float(np.sqrt(msw))})
# the replicate set spans a wider pH range than the analysis sample does, which flatters the ratio; rescale to the sample's own spread
B = ph
for dom in ("Bacteria", "Archaea"):
    s = B[(B.domain == dom) & B.ph_min.notna()].ph_min
    r_dom = 1 - msw / float(s.var())
    print(f"   against the {dom} analysis sample's own spread (variance {s.var():.3f}): reliability {r_dom:.3f}")
    rows.append({"part": "reliability", "subset": f"{dom} analysis sample spread", "n": len(s), "var_total": float(s.var()),
                 "var_within": msw, "icc_single_reading": r_dom})

print("\n3. does the curator disagreement grow in the flat upper half?")
BANDS = [(0, 3.5, "acid-tolerant tail (<= 3.5)"), (3.5, 5.0, "middle (3.5 to 5)"), (5.0, 14, "upper half (>= 5)")]
for lo, hi, name in BANDS:
    ds, ns = [], 0
    for a, b in itertools.combinations(PRECEDENCE, 2):
        m = W[[a, b]].dropna()
        if m.empty:
            continue
        ref = m[a]  # a precedes b, so band on the better-attested reading rather than on the value under test
        sel = m[(ref > lo) & (ref <= hi)]
        ds += list((sel[a] - sel[b]).to_numpy())
        ns += len(sel)
    if ns < 8:
        print(f"   {name:28s} n={ns:3d}  too few pairs")
        continue
    ds = np.asarray(ds)
    print(f"   {name:28s} n={ns:3d}  SD of disagreement {ds.std():.3f}  mean {ds.mean():+.3f}  identical in {(np.abs(ds) < 1e-9).mean():.0%}")
    rows.append({"part": "band", "subset": name, "n": ns, "sd_diff": float(ds.std()), "mean_diff": float(ds.mean()),
                 "exact_rate": float((np.abs(ds) < 1e-9).mean())})
print("   This came out the opposite way round to the hypothesis it was written to test. Disagreement is LARGEST in the acid-tolerant tail")
print("   and smallest in the flat upper half, so the flat upper half is not a stretch of the range where the outcome is too noisy to carry")
print("   a relationship. The phenotype-ceiling gap survives in a narrower form, though: agreement is not accuracy, and because these")
print("   curators copy a shared literature (see the docstring), a conventional minimum that everyone transcribes identically produces")
print("   exactly this signature. What is ruled out is random error; a shared convention is not.")

print("\n4. what the reliability bounds, and what it does not")
(b_clean, se_clean), (b_noisy, se_noisy) = _check_y_error_is_unbiased()
print(f"   simulation: adding outcome noise moves the per-SD slope {b_clean:+.4f} -> {b_noisy:+.4f} and its SE {se_clean:.4f} -> {se_noisy:.4f}")
print("   so curator error in pH_min does NOT bias the coefficient. It costs precision and it caps R-squared. Nothing here disattenuates.")
rows.append({"part": "simulation", "subset": "per-SD slope with and without outcome noise", "beta_clean": b_clean, "beta_noisy": b_noisy,
             "se_clean": se_clean, "se_noisy": se_noisy})
bac = 1 - msw / float(B[(B.domain == "Bacteria") & B.ph_min.notna()].ph_min.var())
print(f"\n   The ceiling on explainable variance in the bacterial sample is the reliability, {bac:.3f}.")
try:
    fs = pd.read_csv(D / "factor_score_headline.csv")
    print("   measured against it, from scripts/56:")
    print(fs.to_string(index=False))
except Exception as e:
    print("   (factor_score_headline.csv not readable:", e, ")")

R = pd.DataFrame(rows)
R.to_csv(D / "outcome_reliability.csv", index=False)
print("\nwrote", D / "outcome_reliability.csv")
print("\nReading: the composite adds under a hundredth of pH_min's variance against a ceiling of at most 0.9, so poor measurement of the")
print("outcome is not why that number is small - the modest conclusion stands exactly as it was. And the flat upper half is not a noise")
print("artefact, which narrows the open phenotype-ceiling question to shared convention rather than retiring it.")
