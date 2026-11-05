"""What the deterministic panel does to the PUBLISHED headline, not to a simplified stand-in -> DATA_DIR/determinism_impact.csv
scripts/79 measured correction 17's impact on a simplified common model (+0.1762 against +0.1772 per SD) because at that point the only comparison available was
rebuild-against-original on a model both builds could be fitted with cheaply. That is the weaker measurement, and quoting it as the impact on the result would be the
same mistake as correction 10: reassurance from a figure that is not the figure anyone relies on. This runs the authoritative table (scripts/57) against the DETERMINISTIC
panel and diffs it row by row against the committed numbers, so the question "does the published headline survive the fix?" is answered on the published model.
THE ANSWER IS YES, AND THE SIZE OF THE MOVE IS THE POINT. The headline goes from +0.1848 to +0.1858 per SD, a shift of 0.0010 against a standard error of about 0.029.
Every one of the eleven model rows keeps its sign, its significance verdict and its ordering; the largest proportional move is in the two controls, which are the rows
with the smallest coefficients and are null before and after (transmembrane -0.0430 -> -0.0460, p 0.073 -> 0.055, still not significant at 0.05; exposed cytoplasmic
+0.0212 -> +0.0191, p 0.33 -> 0.38). Two derived shares move by one point (optimum 57% -> 56%, and its complement) and two ratios in the last decimal. The measured move,
0.035 of a standard error, lands on top of scripts/79's independently-computed 0.036 from a different model, which is the agreement one wants between a cheap
stand-in and the real thing. A built-in placebo guards the comparison: the three hand-picked-panel rows read a file this script does not swap and must move by
exactly zero, and the script fails if they do not.
WHY THE COMMITTED PANEL IS STILL CANONICAL, STATED SO IT IS A DECISION RATHER THAN AN OVERSIGHT. Promoting the deterministic build would require re-running every
downstream script that reads the panel, because a document quoting a new headline beside result files built from the old panel is worse than either panel alone - that
mixture is unauditable, and scripts/69 would not catch it since both figures exist somewhere. The migration is therefore all-or-nothing and is recorded as open work in
docs/SYNTHESIS.md section 8 rather than done halfway. What this script establishes is the thing that decision needs: the cost of NOT migrating is bounded at 0.001 on
the headline and no qualitative claim, which is why it can wait, and the cost is MEASURED rather than assumed.
Usage: python3 scripts/80_determinism_impact.py   (env DATA_DIR=results/broad2; needs structure_panel_wide_det.csv)"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import pandas as pd

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
CANON, DET = D / "structure_panel_wide.csv", D / "structure_panel_wide_det.csv"
NUM, HEAD = ROOT / "docs/NUMBERS.md", D / "headline_numbers.csv"
for f in (CANON, DET):
    if not f.exists():
        sys.exit(f"missing {f}")

print("1. refitting the authoritative table on the deterministic panel")
before = pd.read_csv(HEAD)
bak = {}
for f in (CANON, NUM, HEAD):                      # restored in the finally block whatever happens
    bak[f] = f.read_bytes()
try:
    shutil.copy(DET, CANON)
    r = subprocess.run([sys.executable, str(ROOT / "scripts/57_headline_table.py")], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"scripts/57 failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    after = pd.read_csv(HEAD)
finally:
    for f, b in bak.items():                      # the committed panel and table are put back before anything else
        f.write_bytes(b)
    print("   committed panel, docs/NUMBERS.md and headline_numbers.csv restored")

print("\n2. every model row, committed against deterministic")
# NAMED, not guessed. The first draft picked the key and coefficient columns with next(...) heuristics, matched
# `n` as the key and `sd_of_predictor` as the coefficient, and reported 18 rows "changing significance" with every
# figure nan - while still printing its closing paragraph asserting that the headline survives. A script that
# states a conclusion its own output contradicts is the exact failure this project keeps correcting, so the columns
# are named and the assertions below refuse to let it happen again.
KEY, BCOL, PCOL, SECOL = "figure", "beta_per_sd", "p", "se_per_sd"
for c in (KEY, BCOL, PCOL, SECOL):
    assert c in before.columns and c in after.columns, f"headline_numbers.csv has no {c} column"
# headline_numbers.csv holds the eleven MODEL rows and nine DERIVED-quantity rows in one file; both carry a
# `figure` label and only the model rows carry a coefficient, so filtering on the label alone let nine empty
# rows into the merge and tripped the assertion below. That is the assertion doing its job.
before = before[before[KEY].notna() & before[BCOL].notna()]
after = after[after[KEY].notna() & after[BCOL].notna()]
m = before[[KEY, BCOL, PCOL]].merge(after[[KEY, BCOL, PCOL]], on=KEY, suffixes=("_old", "_new"))
assert len(m) >= 10, f"expected the eleven model rows; matched {len(m)}"
assert m[[f"{BCOL}_old", f"{BCOL}_new"]].notna().all().all(), "a coefficient came back missing"
key, bcol, pcol = KEY, BCOL, PCOL
print(f"   {'row':62s} {'committed':>10s} {'determ.':>10s} {'move':>8s}  sig same")
rows = []
for x in m.itertuples():
    bo, bn = getattr(x, f"{bcol}_old"), getattr(x, f"{bcol}_new")
    po, pn = getattr(x, f"{pcol}_old"), getattr(x, f"{pcol}_new")
    same = (po < 0.05) == (pn < 0.05) and (bo > 0) == (bn > 0)
    print(f"   {str(getattr(x, key))[:62]:62s} {bo:>+10.4f} {bn:>+10.4f} {bn - bo:>+8.4f}  {'yes' if same else 'NO'}")
    rows.append({"part": "model row", "row": getattr(x, key), "beta_committed": bo, "beta_deterministic": bn,
                 "move": bn - bo, "p_committed": po, "p_deterministic": pn, "sign_and_significance_same": same})
flipped = [r for r in rows if not r["sign_and_significance_same"]]
print(f"\n   rows changing sign or crossing p=0.05: {len(flipped)}"
      + ("".join(f"\n      {r['row']}" for r in flipped) if flipped else " - none"))
print(f"   largest absolute move: {max(abs(r['move']) for r in rows):.4f}")
# A PLACEBO THAT COSTS NOTHING, and the comparison is worth less without it. The three hand-picked-panel rows read
# structure_panel.csv, which this script does not touch, so they must move by EXACTLY zero. If they moved, the swap
# would be reaching something it should not and the whole table's movement would be uninterpretable.
placebo = [r for r in rows if "hand-picked" in str(r["row"])]
worst = max((abs(r["move"]) for r in placebo), default=None)
print(f"   placebo: {len(placebo)} hand-picked-panel rows read a panel this script does not swap; largest move "
      f"{worst:.4f}" + (" - exactly zero, as required" if worst == 0 else " - NONZERO, the swap is reaching too far"))
assert placebo and worst == 0, "a row that cannot depend on the swapped panel moved anyway"
rows.append({"part": "placebo", "row": "hand-picked panel rows, which read an unswapped file",
             "n_rows": len(placebo), "largest_move": worst})

print("\n3. the headline specifically, against its own standard error")
h = [r for r in rows if "rule-selected" in str(r["row"])]
if h:
    r0 = h[0]
    se = float(before.loc[before[key] == r0["row"], SECOL].iloc[0])
    print(f"   +{r0['beta_committed']:.4f} -> +{r0['beta_deterministic']:.4f}, a move of {abs(r0['move']):.4f}"
          + (f" = {abs(r0['move']) / se:.3f} of one standard error ({se:.4f})" if se == se else ""))
    rows.append({"part": "headline", "row": r0["row"], "move": r0["move"], "se": se,
                 "move_in_standard_errors": abs(r0["move"]) / se if se == se else None})
pd.DataFrame(rows).to_csv(D / "determinism_impact.csv", index=False)
print("\nwrote", D / "determinism_impact.csv")
# The closing paragraph below states a conclusion. These assertions are what entitle it to: if a row ever does flip,
# or the headline moves by more than a quarter of its standard error, the script fails instead of reassuring.
assert not flipped, f"{len(flipped)} rows changed sign or crossed p=0.05 - the reading below is no longer true"
assert h and abs(h[0]["move"]) < 0.25 * se, "the headline moved more than a quarter of its standard error"
print("\nReading. The published headline survives the determinism fix: +0.1848 becomes +0.1858, well under a tenth of")
print("its own standard error, and no row changes sign or significance. The committed panel stays canonical because")
print("migrating to the deterministic one means re-running everything that reads it - a document quoting the new")
print("headline next to result files built from the old panel would be a mixture neither panel justifies, and")
print("scripts/69 could not catch it because both numbers exist somewhere. That migration is open work, and this")
print("script is what makes deferring it a measured decision instead of an assumption.")
