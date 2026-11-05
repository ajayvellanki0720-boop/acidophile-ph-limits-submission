"""Does every number quoted in this project's documents still exist in its results? -> DATA_DIR/document_consistency.csv
Correction 12 was a figure of +0.042 per SD that one section of docs/SYNTHESIS.md recorded as withdrawn while another went on quoting it as current. It
survived for weeks because nothing ever read the two sections against each other, and because docs/NUMBERS.md - which exists to prevent exactly this -
held both of the estimates involved but not the ratio between them, and the ratio was the claim that was wrong. The fix applied there was to compute
every derived quantity inside the numbers table. This is the general remedy: check mechanically whether each coefficient a document asserts can still be
found in the result file that the same claim cites.
HOW IT GETS ITS POWER, WHICH IS THE WHOLE DESIGN. The obvious version of this check is worthless. Matched against every numeric value in every result
file - a quarter of a million of them - a MADE-UP coefficient matches 100% of the time at two decimal places, 84% at three and 70% at four, and most
figures in these documents are quoted to two or three. A check that cannot fail is not a check. What rescues it is a convention the project already
follows for unrelated reasons: every script declares the file it writes in the first line of its docstring, and the documents cite scripts inline beside
the claims they support. So each figure is matched against the file ITS OWN CLAIM cites, which shrinks the candidate pool from 250,000 values to a few
dozen and drops the chance-match rate to roughly 22%, 7% and 2% at two, three and four decimal places. Both rates are measured on every run, so a clean
result is never read without also seeing what it is worth.
Three design decisions follow from that, and each is reported rather than silent.
  - The unit of a claim is the BULLET, not the blank-line paragraph. Section 6 of the synthesis is a bullet list with no blank lines in it, so splitting
    on blank lines collapsed nine bullets citing nine different scripts into one unit and matched each bullet's figures against the others' files.
  - THE CHECK'S POWER IS REPORTED AS A COUNT, not only as a rate, because a rate is easy to read as reassurance and its
    weight depends entirely on how precisely the documents happen to quote. Part 3 multiplies the per-precision chance-match
    rate by the number of figures actually quoted at that precision, giving the expected number of stale figures a clean run
    would pass. It stands at about 45 of 517, i.e. 9%, which is the number to quote when calling a run clean, having
    been 42 of 492 before 1-dp figures were calibrated and 44 before the five-decimal extension brought 25 further
    figures into the check. It is a SNAPSHOT and moves whenever a document gains or loses figures, so part 3
    reprints it on every run and the run's number, not this line, is the one to use. The single largest weakness it
    names is small and fixable: three figures still quoted to one decimal place account for 1.5 of the 45.
  - AND IT RECORDS A TRADE-OFF RATHER THAN HIDING IT. Declaring 11 previously-undeclared outputs (scripts/77 check C3) was
    the right fix - an undeclared output makes any figure in it unreachable, which is a guaranteed false miss - but it also
    enlarged the candidate pools and so WEAKENED this check, from 38 expected-missed figures to 42. Correctness of the
    declarations and power of the check pull against each other, both are measured, and neither is traded silently.
  - A claim citing more than MAX_CITES scripts has no constraining pool, since the union of many files approaches the global pool, so it is counted as
    unsourced rather than passed on a near-certain match. The README contains one historical line citing 34 scripts.
  - Deliberately quoted historical figures are not errors. The synthesis's corrections section exists to quote withdrawn numbers and is excluded by
    heading; withdrawn figures quoted elsewhere are listed in KNOWN_HISTORICAL with what retired them.
WHAT A CLEAN RUN MEANS IS NARROW. Every coefficient asserted as current traces to a stored result. It does not mean the figure is the RIGHT one for the
claim it is attached to - no script judges that - and it cannot see a claim that is wrong without being numeric.
On its first full run this found four things: a stale disattenuation range in the README that correction 7 had superseded, where the synthesis was fixed
at the time and the README was not; a missing citation on the one-sidedness paragraph; two correlations that scripts/54 printed and never stored, so the
documents had been quoting them from a terminal; and a withdrawn figure that needed annotating rather than fixing.
Usage: python3 scripts/69_document_consistency.py   (env DATA_DIR=results/broad2)"""
import glob
import os
import re
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
# docs/DATABASES.md joined this list on 2026-11-01, when scripts/94 put measured coefficients into it for the first
# time. Until then it held no figures from this project's own fits, only counts describing other people's databases,
# and leaving it out was harmless. The moment a coefficient was written there it became exactly what correction 12
# was: a number in a document nothing checks, free to outlive the analysis that produced it.
DOCS = ["docs/SYNTHESIS.md", "docs/CHARGE_REVERSAL_DESIGN.md", "docs/PREREGISTRATION.md", "README.md", "docs/DATABASES.md",
        "docs/MANUSCRIPT.md", "docs/EXPERIMENTAL_EVIDENCE.md"]   # joined 2026-11-01; it is the document that will leave the project, so it is the one that least tolerates a stale figure
SKIPPED = {"docs/NUMBERS.md": "generated by scripts/57, so it cannot disagree with the results",
           "docs/LIVING_RESEARCH_DOCUMENT.md": "a dated record that cites no scripts; its figures were correct on their date"}
EXCLUDE_HEADINGS = ("## 5. Corrections to my own analyses",)
MAX_CITES = 4            # above this the union of cited files approaches the global pool and the check loses its power
# FIVE DECIMAL PLACES, not four. The cap was 4, and the trailing (?![\d]) means a figure quoted to FIVE places
# did not match at all - it was not checked and was not counted as unchecked either, so it was invisible to this
# script in both directions. scripts/82, 83 and 84 quote per-point coefficients to five places, and adding a
# README entry with fifteen such figures raised this script's README count by one. The advice part 3 prints -
# quote to four places rather than two, because the same check is nine times stronger there - therefore had a
# cliff immediately past its own recommendation, where quoting MORE precisely removed a figure from the check
# entirely. Five places is calibrated with the rest below.
# SIX decimals since 2026-11-01, five before that and four before the correction-14 work. The cap is a BLIND SPOT, not
# a formatting preference: a figure quoted to more places than this matches nothing, is counted nowhere, and reports as
# neither checked nor untraceable - it is simply invisible. That is the fifth time this family has appeared, and it
# surfaced because scripts/95 reproduces the headline to +0.184838 against +0.184837, an agreement that cannot be stated
# in five places at all. Raising the cap is the fix; noticing that a quiet non-match is indistinguishable from a pass
# is the lesson, and it is why 2d counts figures covered NOWHERE rather than trusting the per-document totals.
COEF = re.compile(r"(?<![\d.])([+−-])(0\.\d{1,6})(?![\d])")
NSIZE = re.compile(r"\bn\s*=\s*(\d[\d,]{1,9})\b")
# sample sizes quoted deliberately although they are not the cited file's own, with why
KNOWN_HISTORICAL_N = {2280: "the hand-picked panel's n, quoted in the README's 2026-10-31 correction note as the wrong value"}
# figures retired by a correction but still quoted outside the corrections section, with what retired them
KNOWN_HISTORICAL = {"+0.042": "correction 6", "+0.07": "correction 6", "+0.071": "correction 4", "+0.11": "correction 9",
                    "+0.109": "correction 9", "+0.217": "correction 10", "+0.52": "the failed method in scripts/65",
                    "+0.158": "correction 7", "+0.192": "correction 7",
                    "+0.225": "an earlier panel's figure, retained deliberately",
                    "+0.101": "an earlier panel's figure, retained deliberately",
                    "+0.3137": "the deliberately bogus value used to mutation-test this very script",
                    # the withdrawn histidine lead's pooled figure, quoted in Section 8 because the pKa test's
                    # predicted direction turns on it. scripts/28 printed it and never stored it, and re-running a
                    # retired analysis to give it a row is not worth it - the lead itself is withdrawn (Section 4).
                    "-0.55": "the withdrawn histidine lead (Section 4); printed by scripts/28, never stored",
                    # A NEW CLASS, and it is worth distinguishing from the rest of this table. Every entry above is
                    # a figure that WAS a result and was later retired. This one never was: scripts/85's first
                    # version estimated Pagel's lambda on a model including the charged-residue controls, giving
                    # -0.13717 where scripts/82 gives -0.13168, and an assertion rejected it before the script ever
                    # wrote a file. The documents quote it to record the error, which is the project's practice, and
                    # that practice produces numbers no result file can contain BY CONSTRUCTION - the script was
                    # fixed before it produced output. This check was right to flag it and will be right to flag
                    # the next one; the annotation is the correct resolution, not a weakening.
                    "-0.13717": "scripts/85's rejected first draft (wrong lambda model); never written to a file",
                    # Same class as the entry above, and the second instance in a day, which is what makes it a
                    # class rather than a one-off: scripts/87's docstring first quoted -0.463 as the correlation
                    # between the exposed acidic and histidine shares. That is the figure on the RAW panel, before
                    # the bacteria-with-a-tip-and-a-pH-minimum filter; the script's own part 1 prints -0.549 on the
                    # sample it fits. The README quotes the wrong value to record the error, so no result file can
                    # contain it. Documenting a rejected draft and keeping every quoted figure traceable pull
                    # against each other, and this table is where that tension is paid rather than hidden.
                    "-0.463": "scripts/87's withdrawn docstring figure (raw panel, pre-filter); never a result",
                    # docs/CHARGE_REVERSAL_DESIGN.md quotes this as the "within-family +0.312" whose apparent
                    # conflict with the between-species figure the same sentence calls an error of mine: a per-unit
                    # slope set against a per-standard-deviation one, which is correction 2. It became checkable
                    # only when that claim's bare `39a` citation was given its scripts/ prefix on 2026-11-01, and
                    # it is annotated rather than removed because the sentence exists to record the error.
                    "+0.312": "correction 2's per-unit/per-SD conflation, quoted to record it",
                    # A THIRD REASON a figure can be unreachable, distinct from the two above: it was a result,
                    # has not been withdrawn, and simply belongs to a run this DATA_DIR does not hold. NuoL's
                    # +0.076 is a per-reference coefficient from the original 7-reference panel; scripts/54's
                    # stored output is the 32-reference run. Its companion +0.088 DOES match that file, which is
                    # worth noticing rather than being pleased about - at three decimals the chance-match rate is
                    # 7.1%, so a figure from the wrong panel passing is exactly what that rate describes. The
                    # document now says so beside both numbers.
                    "+0.076": "a 7-reference-panel per-reference coefficient; that run is not stored in this DATA_DIR",
                    # A FOURTH REASON, and the mildest: a figure quoted only to record that it USED to be written
                    # that way. docs/CHARGE_REVERSAL_DESIGN.md says the slope range's endpoints "were previously
                    # rounded to +0.10 and +0.19", and the point of the sentence is that a rounded endpoint
                    # matches no stored value and so cannot be checked - which is why the exact figures replaced
                    # it. The rounded pair is kept in the prose to record the change and annotated here.
                    "+0.10": "a rounded range endpoint, superseded in place by +0.101; kept to record the change",
                    "+0.19": "a rounded range endpoint, superseded in place by +0.185; kept to record the change"}


def rnd(x, dp):
    """Round HALF-UP, the way a document's f-string does, and from full precision.

    THIS FIXES A DEFECT THAT EXTENDING THE CHECK TO FIVE DECIMALS EXPOSED. The pools used to be stored at six
    decimals and matched with Python's round(), which is half-to-EVEN. A coefficient printed by `f"{v:.5f}"` as
    +0.11187 came from 0.1118653906...; stored at six decimals that is exactly 0.111865, an exact half at five
    decimals, and half-to-even sent it to 0.11186. The figure in the document was right, the value in the result
    file was right, and the check called it untraceable - a FALSE MISS manufactured by rounding twice. It could
    only bite at five decimals, because that is the only precision one digit from the six-decimal store, which is
    why four years of this check at two to four decimals never showed it.
    """
    f = 10.0 ** dp
    return float(np.floor(np.abs(x) * f + 0.5) / f) * (1 if x >= 0 else -1)


# A REGRESSION GUARD THAT RUNS EVERY TIME, on the exact value that was wrong. scripts/84 reports the acidic
# exposed-periplasmic coefficient under its scoring-quality control as +0.11187, from 0.1118653906829323. Under the
# old two-stage rounding this check called that figure untraceable while both the document and the result file were
# correct. The assertion is on the VALUE rather than on rnd's general behaviour because a general property is easy
# to satisfy accidentally, and because this is the case that actually failed.
assert rnd(0.1118653906829323, 5) == 0.11187, "half-up rounding is broken; five-decimal figures will false-miss"
assert rnd(np.round(0.1118653906829323, 6), 5) == 0.11187, "an exact half at 5 dp must round up, not to even"
assert rnd(-0.1118653906829323, 5) == -0.11187, "negative figures must round away from zero, as a document prints"
assert rnd(0.1848374999999, 6) == 0.184837 and rnd(0.1848375000001, 6) == 0.184838, "six-decimal rounding must be half-up on the figure scripts/95 reproduces"


_ROUNDED = {}


def rounded(pool, dp):
    """The pool rounded to `dp`, built once per (pool, precision).

    PURELY A SPEEDUP, and it was needed. Both the per-figure check and part 3's calibration rebuilt
    `{rnd(x, dp) for x in pool}` inside their innermost loop - the calibration does 400 draws against each of 131
    pools at each of five precisions, so it rebuilt the same rounded set 400 times over. With the global pool at
    354,000 values a full run had grown past fifteen minutes, which is long enough that it stops being run after
    every edit, and a check that is not run is not a check. Results are identical by construction.
    """
    k = (id(pool), dp)
    if k not in _ROUNDED:
        _ROUNDED[k] = {rnd(x, dp) for x in pool}
    return _ROUNDED[k]


def numeric_values(path):
    vals = set()
    try:
        d = pd.read_csv(path)
    except Exception:
        return vals
    for c in d.columns:
        v = pd.to_numeric(d[c], errors="coerce").dropna()
        # NOT rounded to 6 here: see rnd() above. An intermediate rounding is what created the
        # exact-half case that made a correct five-decimal figure untraceable. 12 places only
        # normalises float representation noise and is far below any precision these documents quote.
        vals |= set(np.round(v.to_numpy(float), 12))
    return vals


def int_values(path):
    """Whole numbers appearing anywhere in a result file - the pool a quoted sample size is checked against.
    Coefficients and sample sizes need separate pools because COEF only matches signed sub-unit decimals, so an
    'n = 2,344' was invisible to this script until 2026-10-31 and a wrong one duly reached the README."""
    out = set()
    try:
        d = pd.read_csv(path)
    except Exception:
        return out
    for c in d.columns:
        v = pd.to_numeric(d[c], errors="coerce").dropna()
        out |= {int(x) for x in v.to_numpy(float) if float(x).is_integer() and 1 < abs(x) < 1e7}
    return out


def claims_of(doc_text):
    """Split into claim units. Bullets, numbered items and TABLE ROWS are each their own unit rather than being
    merged into the blank-line block containing them. Section 6 of the synthesis is a bullet list with no blank
    lines, and the evidence table is a markdown table whose rows cite different scripts - treating either as one
    unit matches every figure in it against every cited file's values, which is how the check loses its power.
    A TABLE ROW INHERITS ITS TABLE'S CITATION. Splitting rows apart is what keeps the check specific, but it also
    threw away the citation, which in these documents sits in the paragraph that introduces the table rather than in
    each row. That left every row of every table unsourced - 30-odd figures in the synthesis alone, checked against
    nothing. A row now carries the citations of the nearest preceding non-table block, unless it names a script of
    its own, so each row is still matched against one file and that file is its own table's."""
    out = []
    inherited = ""
    for blk in re.split(r"\n\s*\n", doc_text):
        if not blk.strip():
            continue
        is_table = all(ln.strip().startswith("|") or not ln.strip() for ln in blk.split("\n"))
        if not is_table:
            found = CITE.findall(blk)
            inherited = " ".join(f"scripts/{n}{' ' + v if v else ''}" for n, v in found)
        start = len(out)                       # only THIS block's rows inherit, not every row added so far
        cur = []
        for line in blk.split("\n"):
            starts_claim = re.match(r"\s*(?:[-*]\s|\d+\.\s|\|)", line)
            if starts_claim and cur:
                out.append("\n".join(cur))
                cur = [line]
            else:
                cur.append(line)
        if cur:
            out.append("\n".join(cur))
        if is_table and inherited:
            # attach the table's own citation to each of its rows, for rows that do not cite a script themselves
            for k in range(start, len(out)):
                if out[k].strip().startswith("|") and "scripts/" not in out[k]:
                    out[k] += f"  [inherited {inherited}]"
    return [c for c in out if c.strip()]


# THE META-SCRIPTS. scripts/69 and scripts/77 analyse other scripts' TEXT, so they contain the panel option's
# own syntax (`== "wide"`) without taking the option. The first attempt to separate them tested whether the
# script reads a structural panel - which failed hilariously, because writing that test put the literal
# "structure_panel" into both checkers' source and they matched themselves. A checker cannot reliably
# classify itself by its own string patterns, so the two are named, with the reason stated.
META_SCRIPTS = {"69": "analyses other scripts' text, so it contains the option syntax without taking it",
                "77": "same"}
PANEL_OPTION = re.compile(r"""==\s*["']wide["']|\("orig",\s*"wide"\)|\bTAG\s*,\s*PANEL_CSV|["']wide["']\s*==""")
CITE = re.compile(r"scripts/(\d+[a-z]?)(?:\s+(wide|orig))?")


def cites_of(claim):
    """Script citations as (number, variant) pairs. `scripts/39a wide` selects that script's _wide output and a bare
    `scripts/39a` its default one, so the two panels' figures cannot be substituted for each other."""
    return sorted({(n, "wide" if v == "wide" else "base") for n, v in CITE.findall(claim)})


def label(cites):
    return ";".join(n + (" wide" if v == "wide" else "") for n, v in cites)


# THIS SCRIPT'S OWN OUTPUT IS EXCLUDED FROM EVERY POOL, and leaving it in was a circularity, not an
# inefficiency. document_consistency.csv records each figure it checked in a `quoted` column - 617 numbers, every
# coefficient these documents assert - so pooling it means a figure is matched against the check's own record of
# that same figure. Any claim citing scripts/69 had its numbers verifying themselves, and several entries written
# on 2026-11-01 cite it. It was found by computing a PER-FIGURE chance-match rate (part 2e): the weakest figures
# in the corpus all cited scripts/69, because its pool is 2,020 values and most of them are the documents'. The
# audit's output is excluded for the same reason - it quotes figures out of the scripts it audits.
SELF_OUTPUT = ("document_consistency.csv", "audit.csv")
print("1. the value pools")
every, nfiles = set(), 0
for f in sorted(glob.glob(str(ROOT / "results/**/*.csv"), recursive=True)):
    if Path(f).name in SELF_OUTPUT:
        continue
    v = numeric_values(f)
    if v:
        every |= v
        nfiles += 1
DPS = tuple(range(1, int(COEF.pattern.split("{1,")[1].split("}")[0]) + 1))   # derived from COEF, so raising the cap cannot leave a precision uncalibrated
GLOBAL = {dp: {rnd(x, dp) for x in every} for dp in DPS}
print(f"   {len(every)} distinct values across {nfiles} result files - the pool that makes the naive version of this worthless")
declared, file_pool, int_pool = {}, {}, {}
for f in sorted(glob.glob(str(ROOT / "scripts/*.py"))):
    head = Path(f).read_text()[:1400]
    num = re.match(r"(\d+[a-z]?)_", Path(f).name)
    # a script may declare SEVERAL outputs on one line ("-> docs/NUMBERS.md and DATA_DIR/headline_numbers.csv"), and
    # taking only the first match silently gave scripts/57 no pool at all, which this check then reported as a miss
    outs = re.findall(r"(?:DATA_DIR|docs)/([A-Za-z0-9_./]+\.csv)", head.split("Usage:")[0]) if num else []
    # AND A LITERAL results/ PATH COUNTS TOO, which it did not until 2026-11-01 - a fifth mechanical blind spot,
    # found while trying to cite the archaea bullet in docs/SYNTHESIS.md. Five scripts declare their output as
    # `results/<dataset>/x.csv` rather than `DATA_DIR/x.csv`, because they write to a fixed dataset instead of a
    # selectable one: 07, 23, 35c, 35d and 43. Every one of those files exists on disk and every one was invisible
    # here, so a claim citing those scripts had no pool at all and its figures were reported unsourced - which is
    # how ten archaeal figures sat unchecked. These resolve from the repository ROOT, not from DATA_DIR, because
    # the dataset is part of what the script declares.
    lits = re.findall(r"\bresults/([A-Za-z0-9_./]+\.csv)", head.split("Usage:")[0]) if num else []
    if num and (outs or lits):
        declared[num.group(1)] = ";".join(outs + [f"results/{o}" for o in lits])
        # PANEL VARIANTS GET SEPARATE POOLS, and this is not a refinement - without it the check cannot do its job.
        # A script taking orig/wide declares BOTH outputs on its first line, so the union of its pool holds the
        # hand-picked panel's figures and the rule-selected panel's alike, and a claim citing it passes on either.
        # That is how a wrong sample size reached the README on 2026-10-31: n=2,280 is the hand-picked panel's and
        # the claim was about the wide one, and scripts/39a's union pool contains both. The documents already
        # write the variant into the citation ("scripts/39a wide"), so the pool is split the same way and a bare
        # citation means the default run - which also stops a bare citation passing on a wide-panel value.
        # ...BUT ONLY FOR SCRIPTS THAT ACTUALLY TAKE THE OPTION. The split guards against a bare citation passing
        # on a wide-panel value, which can only happen where one script produces both panels. scripts/78 takes no
        # panel argument and runs only on the rule-selected references, so its single output is NAMED
        # distribution_wide.csv - and the naive split filed its whole pool under "wide", leaving the "base" pool
        # empty and every bare `scripts/78` citation an automatic miss. That cost two real figures their check
        # and would have been read as the document being wrong, when the document was right and this was. A
        # script with no option has no ambiguity to guard against, so its outputs go into both pools.
        # Detected on the WHOLE file, not on `head`. The option is selected in the CODE and `head` is only the
        # first 1400 characters, i.e. the docstring - so detecting it there found none of the four scripts that
        # really take it (36c, 39a, 39b, 75) and would have merged both panels' pools for exactly the scripts the
        # split exists for, destroying the guard that caught the n=2,280 error. Only 45 and 78 are genuinely
        # option-less with a _wide-named output.
        # The pattern must match how each script ACTUALLY branches, not one spelling of it. The first version
        # looked for `sys.argv[1] ... wide` on a single line and missed scripts/39b, which reads the argument on
        # one line and compares on another (`if WHICH == "wide":`), so 39b's two panels would have been merged.
        # Matching the COMPARISON catches every form. scripts/77's check C4 carried the same narrow pattern and
        # was mis-detecting 39b too; both now use this one.
        # ...and not one of the meta-scripts, which contain the option's syntax without taking it (see
        # META_SCRIPTS above). Merging their pools would do no harm here, but the same test serves scripts/77's
        # check C4, where it does.
        takes_panel = bool(PANEL_OPTION.search(Path(f).read_text())) and num.group(1) not in META_SCRIPTS
        split = (("wide", [o for o in outs if o.endswith("_wide.csv")]),
                 ("base", [o for o in outs if not o.endswith("_wide.csv")])) if takes_panel else \
                (("wide", outs), ("base", outs))
        for key, sel in split:
            v, iv = set(), set()
            for o in sel:
                if Path(o).name in SELF_OUTPUT:         # circular; see the comment at SELF_OUTPUT
                    continue
                v |= numeric_values(D / o)
                iv |= int_values(D / o)
            for o in lits:                       # resolved from the root; see the comment at `lits` above
                v |= numeric_values(ROOT / "results" / o)
                iv |= int_values(ROOT / "results" / o)
            if v:
                file_pool[(num.group(1), key)] = v
            if iv:
                int_pool[(num.group(1), key)] = iv
print(f"   {len(declared)} scripts declare an output file; {len(file_pool)} script-variant pools exist, median "
      f"{int(np.median([len(v) for v in file_pool.values()]))} values each")
print(f"   {sum(1 for k in file_pool if k[1] == 'wide')} of those pools are a panel variant kept separate from the "
      f"same script's default output")

rows, all_missing, size_checks = [], [], []
for name, why in SKIPPED.items():
    print(f"   skipping {name}: {why}")
    rows.append({"document": name, "quoted": "SKIPPED", "note": why})

print("\n2. each coefficient against the file its own claim cites")
figrates = []

for doc in DOCS:
    text = (ROOT / doc).read_text()
    for h in EXCLUDE_HEADINGS:
        i = text.find(h)
        if i < 0:
            continue
        j = text.find("\n## ", i + len(h))
        text = text[:i] + text[len(text) if j < 0 else j:]
    checked = unsourced = over = 0
    missing = []
    nchecked = nunsourced = 0
    for claim in claims_of(text):
        # SAMPLE SIZES, checked the same way and reported separately. COEF only matches signed sub-unit decimals,
        # so until 2026-10-31 an "n = 2,344" was invisible here; a claim about the wide panel quoting the hand-picked
        # panel's n=2,280 passed every check this script ran.
        sizes = sorted({int(m.replace(",", "")) for m in NSIZE.findall(claim)})
        if sizes:
            ncites = cites_of(claim)
            npools = [int_pool[c] for c in ncites if c in int_pool]
            if len(npools) > MAX_CITES:
                npools = []
            for sz in sizes:
                nnote = KNOWN_HISTORICAL_N.get(sz, "")
                if not npools:
                    nunsourced += 1
                    rows.append({"document": doc, "quoted": f"n={sz}", "kind": "sample size",
                                 "cites": label(ncites), "scope": "no cited script",
                                 "found": any(sz in p for p in int_pool.values()), "note": nnote})
                    continue
                nchecked += 1
                size_checks.append((sz, npools))
                hit = any(sz in p for p in npools)
                rows.append({"document": doc, "quoted": f"n={sz}", "kind": "sample size",
                             "cites": label(ncites), "scope": "cited script's own file", "found": hit, "note": nnote})
                if not hit and not nnote:
                    missing.append((f"n={sz}", label(ncites), re.sub(r"\s+", " ", claim.strip())[:100]))
        figs = sorted({("-" if sg in "-−" else "+", nm) for sg, nm in COEF.findall(claim)})
        if not figs:
            continue
        cites = cites_of(claim)
        pools = [file_pool[c] for c in cites if c in file_pool]
        too_many = len(pools) > MAX_CITES
        if too_many:
            over += len(figs)
            pools = []
        for sign, num in figs:
            dp = len(num.split(".")[1])
            val = float(num) * (1 if sign == "+" else -1)
            key = f"{sign}{num}"
            note = KNOWN_HISTORICAL.get(key, "")
            if not pools:
                unsourced += 1
                rows.append({"document": doc, "quoted": key, "decimals": dp, "cites": label(cites),
                             "scope": "too many citations to constrain" if too_many else "no cited script",
                             "found": rnd(val, dp) in GLOBAL[dp], "note": note})
                continue
            checked += 1
            hit = any(rnd(val, dp) in rounded(pl, dp) for pl in pools)
            # THE PER-FIGURE CHANCE-MATCH RATE, recorded for every figure that passes. Part 3's corpus average
            # is the right number for "how much is a clean run worth", and it is the wrong number for "is THIS
            # claim verified": a figure quoted to two decimals against a 200-value pool is barely constrained,
            # while the same figure to five decimals against twenty values is strong, and the average hides
            # which is which. It is the fraction of representable values in (-0.6, 0.6) that this figure's own
            # pools occupy at its own precision. It is what found the circularity fixed on 2026-11-01: the
            # weakest figures in all four documents cited scripts/69, whose output recorded the documents' own
            # figures, so they were matching themselves.
            per_fig = (len({z for pl in pools for z in rounded(pl, dp) if -0.6 < z < 0.6}) / (1.2 * 10 ** dp)
                       if hit else None)
            if per_fig is not None:
                figrates.append((per_fig, doc, key, label(cites), dp))
            rows.append({"document": doc, "quoted": key, "decimals": dp, "cites": label(cites),
                         "scope": "cited script's own file", "found": hit, "note": note,
                         "chance_match_this_figure": per_fig})
            if not hit and not note:
                missing.append((key, label(cites), re.sub(r"\s+", " ", claim.strip())[:100]))
    print(f"   {doc:36s} {checked:4d} checked, {unsourced:4d} unsourced ({over} of them over the citation limit), "
          f"{len(missing)} untraceable | sample sizes: {nchecked} checked, {nunsourced} unsourced")
    rows.append({"document": doc, "quoted": "DOCUMENT SUMMARY", "n_checked": checked, "n_unsourced": unsourced,
                 "n_over_cited": over, "n_untraceable": len(missing),
                 "n_sizes_checked": nchecked, "n_sizes_unsourced": nunsourced})
    all_missing += [(doc, *m) for m in missing]

print("\n2c. exposure to the MECHANICAL false-miss modes, which this check's calibration does not measure")
# The heading used to say "the two" modes and then listed four, because two more were found after it was written.
# That is check C6's error class - a prose count contradicting what it counts - in the script whose own report
# carries it, so the number is simply not stated rather than stated and left to drift.
# WHY THIS EXISTS. Part 3 measures how often a WRONG figure would pass - the false-negative rate - with some care.
# It says nothing about the opposite error, a RIGHT figure reported untraceable, and on 2026-11-01 this check
# produced three of those in one day of document edits. Each had a mechanical cause and none would have been
# caught by reading: a citation naming a panel-option script without its variant, so the pool was built from a
# base spelling whose file does not exist; a figure quoted in order to RECORD a rejected draft, which no result
# file can contain by construction; and an input a script consumed to compute a ratio without ever writing it to
# its own output. The first two are countable in advance, which is what this part does. A high count here does not
# mean anything is wrong - it means that many figures are one edit away from being unreachable.
missing_decl = []
for num, outs in sorted(declared.items()):
    for o in outs.split(";"):
        if not (D / o).exists():
            missing_decl.append((num, o))
opt_scripts = {n for n in declared
               if (ROOT / "scripts").glob(f"{n}_*.py") and n not in META_SCRIPTS
               and any(PANEL_OPTION.search(f.read_text()) for f in (ROOT / "scripts").glob(f"{n}_*.py"))}
# COUNTED TWICE, because the raw count overstates the risk. A bare citation only matters if a figure is quoted
# beside it: a prose mention of a script with no number attached cannot produce a false miss. The second count is
# the actionable one, and splitting them is the difference between a number that prompts a fix and one that
# prompts 34 edits of which most change nothing.
bare, bare_live = [], []
for doc in DOCS:
    body = (ROOT / doc).read_text()
    for n, v in CITE.findall(body):
        if n in opt_scripts and v != "wide":
            bare.append((doc, n))
    for claim in re.split(r"\n\s*\n|\n(?=[-*|] )", body):
        if not COEF.search(claim):
            continue
        for n, v in CITE.findall(claim):
            if n in opt_scripts and v != "wide":
                bare_live.append((doc, n))
print(f"   {len(missing_decl)} declared outputs do not exist in this DATA_DIR, so every figure in them is")
print("   unreachable: " + (", ".join(f"{n}:{o}" for n, o in missing_decl[:6]) + ("..." if len(missing_decl) > 6 else "")
                            if missing_decl else "(none)"))
# AND NARROWED ONCE MORE, to the set that will actually bite. A bare citation beside a figure is harmless where
# the script's BASE output exists, which is why 30 of them pass today: most of these scripts have been run on both
# panels. The risk is specific to scripts whose base output is missing - there, the pool is empty and the next
# figure quoted from the _wide file fails silently. That intersection is the standing debt, and it is small enough
# to fix by hand, which a count of 34 was not.
# TESTED ON THE POOL, NOT ON THE FILE LIST, because the first version got this wrong and overstated the risk.
# It treated a script as having no base pool if ANY of its non-_wide declared outputs was missing, which flagged
# scripts/36c - that script declares several outputs, one of them absent, and the rest give it a perfectly good
# base pool. The condition that matters is whether the base pool EXISTS at all, which file_pool already answers.
# The wrong version named 13 at-risk citations and the right one names fewer; a count that is too high is not the
# safe direction to be wrong in, because it invites 13 edits of which several would have introduced false misses.
no_base = {n for n in opt_scripts if (n, "base") not in file_pool}
at_risk = sorted({(doc, n) for doc, n in bare_live if n in no_base})
print(f"   {len(bare)} citations name a panel-option script WITHOUT a variant, and {len(bare_live)} of those sit in")
print("   a claim that quotes a coefficient. Most are harmless because the script's base output exists; the ones")
print(f"   that will bite are the {len(at_risk)} where it does not, because there the pool is empty and the next")
print("   figure quoted from the _wide file fails silently:")
print("      " + (", ".join(f"{d.split('/')[-1]}:{n}" for d, n in at_risk) if at_risk else "(none)"))
print("   This was the cause of two of the three false misses on 2026-11-01.")
# A FOURTH MODE, found on 2026-11-01 by reading the design document rather than by any check. Several claims cite
# their script in a BARE form - `39a` rather than `scripts/39a` - which the citation pattern does not match, so
# the claim looked uncited and every figure in it was counted unsourced. docs/CHARGE_REVERSAL_DESIGN.md had 39
# figures checked and 64 unsourced; prefixing two citations took it to 55 checked, and the nine figures that then
# failed were real - two belonged to scripts/57, one was a withdrawn figure needing annotation, and two were
# range endpoints rounded to +0.10 and +0.19 where the stored values are +0.101 and +0.185, so they matched
# nothing. An uncited claim is not a passing claim; it is an unexamined one, and this mode makes a claim look
# uncited while its author believed it was sourced.
KNOWN_SCRIPTS = {re.match(r"(\d+[a-z]?)_", Path(f).name).group(1)
                 for f in glob.glob(str(ROOT / "scripts/*.py"))
                 if re.match(r"(\d+[a-z]?)_", Path(f).name)}
BARE_CITE = re.compile(r"`(\d+[a-z]?)`")
bare_prefix = []
for doc in DOCS:
    for claim in re.split(r"\n\s*\n|\n(?=[-*|] )", (ROOT / doc).read_text()):
        if not COEF.search(claim) or CITE.search(claim):
            continue
        for b in BARE_CITE.findall(claim):
            if b in KNOWN_SCRIPTS:
                bare_prefix.append((doc, b))
print(f"   {len(bare_prefix)} figure-bearing claims cite a script WITHOUT the scripts/ prefix, so the citation is")
print("   invisible here and the claim looks uncited: " + (", ".join(f"{d.split('/')[-1]}:{b}" for d, b in bare_prefix[:8]) if bare_prefix else "(none)"))
rows.append({"document": "ALL", "quoted": "BARE PREFIXLESS CITATIONS", "n_bare_prefixless": len(bare_prefix),
             "where": ";".join(f"{d}:{b}" for d, b in bare_prefix[:30])})
print(f"   {len(KNOWN_HISTORICAL)} figures are exempted by annotation rather than found, which is the third mode")
print("   paid for honestly: a number quoted to record an error cannot be in a result file.")
rows.append({"document": "ALL", "quoted": "FALSE MISS EXPOSURE",
             "n_declared_outputs_missing": len(missing_decl), "n_bare_panel_citations": len(bare), "n_bare_panel_citations_beside_a_figure": len(bare_live), "n_bare_citations_at_risk": len(at_risk),
             "n_annotated_exemptions": len(KNOWN_HISTORICAL)})

print("\n2b. calibration of the sample-size check, against the two ways a wrong n actually arises")
if size_checks:
    rng_n = np.random.default_rng(7)
    quoted = [s for s, _ in size_checks]
    near = swap = 0
    for sz, pls in size_checks:
        # (a) a nearby wrong number: a typo, or an n read off the wrong row of the right file
        off = 0
        while off == 0:
            off = int(rng_n.integers(-300, 301))
        near += any((sz + off) in p for p in pls)
        # (b) ANOTHER n quoted somewhere in these documents, substituted for this one. This is the error that
        # occurred: the wide panel's claim was given the hand-picked panel's n, and both are real sample sizes.
        others = [q for q in quoted if q != sz]
        if others:
            swap += any(int(rng_n.choice(others)) in p for p in pls)
    print(f"   {len(size_checks)} sample sizes checked against a cited file")
    print(f"   a nearby wrong value (+/-1 to 300) still matches {near / len(size_checks):.0%} of the time")
    print(f"   another n quoted in these documents, substituted, matches {swap / len(size_checks):.0%} of the time")
    print("   The second rate is the one that matters, because that substitution is what happened; it is low only")
    print("   because panel variants get separate pools. Without that split the same rate would be near certain for")
    print("   every script that takes an orig/wide argument, since both panels' outputs are declared on one line.")
    rows.append({"document": "CALIBRATION", "quoted": "sample size", "kind": "sample size",
                 "n_checked": len(size_checks), "miss_rate_nearby": near / len(size_checks),
                 "miss_rate_swapped_quoted_n": swap / len(size_checks)})
else:
    print("   no sample size in any document cites a script whose output file exists - nothing to calibrate")

print("\n2d. does a document's uncovered figure matter? only if it is uncovered EVERYWHERE")
# WHY THIS IS NOT PEDANTRY. The README reports 70 unsourced figures, 68 of them over the citation limit, which
# reads as a large hole. But these documents overlap heavily: a figure the README quotes inside a long changelog
# entry citing eight scripts is usually the same figure docs/SYNTHESIS.md quotes beside one or two, where it IS
# checked. A per-document count therefore overstates the gap, and the number that matters is how many distinct
# figures are checked in NO document at all. Reporting only the per-document count would be the error correction
# 10 was about in miniature: a proxy quoted where the quantity itself is available.
checked_ok = {str(r["quoted"]) for r in rows if r.get("found") is True and r.get("scope") == "cited script's own file"}
uncovered = {}
for r in rows:
    if r.get("kind") == "sample size" or r.get("decimals") is None:
        continue
    key = str(r.get("quoted", ""))
    if r.get("scope") in ("too many citations to constrain", "no cited script") and key not in checked_ok:
        uncovered.setdefault(key, set()).add(r["document"])
per_doc = sum(1 for r in rows if r.get("decimals") is not None and r.get("kind") != "sample size"
              and r.get("scope") in ("too many citations to constrain", "no cited script"))
print(f"   {per_doc} figure mentions are unsourced across the four documents, but they are not {per_doc} distinct")
print(f"   problems: {len(uncovered)} DISTINCT figures are checked in no document at all. The rest are quoted")
print("   somewhere that constrains them, usually docs/SYNTHESIS.md beside one or two citations rather than the")
print("   README inside a changelog entry citing eight.")
if uncovered:
    sample = sorted(uncovered)[:10]
    print("   uncovered everywhere: " + ", ".join(sample) + ("..." if len(uncovered) > 10 else ""))
# AND A STRUCTURAL BLIND SPOT WORTH NAMING, found by asking where the uncovered figures live. The corrections
# section of docs/SYNTHESIS.md is excluded from checking BY DESIGN, because it exists to quote withdrawn numbers
# and would otherwise fail on every one. The consequence nobody had drawn: a figure introduced only in a
# correction entry is checked NOWHERE - the section that holds it is exempt, and the README's parallel changelog
# entry cites too many scripts to constrain anything. Corrections are exactly where new figures get introduced,
# so the exemption that makes the check possible also creates its largest hole. The remedy is per-figure and not
# structural: a correction's CURRENT figures should also appear in a checkable section beside a tight citation.
synth = (ROOT / "docs/SYNTHESIS.md").read_text()
corr = synth[synth.index(EXCLUDE_HEADINGS[0]):] if EXCLUDE_HEADINGS[0] in synth else ""
corr = corr[:corr.index("\n## ", 10)] if "\n## " in corr[10:] else corr
only_corr = [k for k in uncovered if k in corr]
print(f"   Of those {len(uncovered)}, {len(only_corr)} appear in the corrections section of docs/SYNTHESIS.md,")
print("   which is exempt from this check BY DESIGN because it exists to quote withdrawn numbers. A figure")
print("   introduced only in a correction entry is therefore checked nowhere: the section holding it is exempt")
print("   and the README's parallel entry cites too many scripts to constrain it. Corrections are where new")
print("   figures get introduced, so the exemption that makes this check possible also creates its largest hole.")
rows.append({"document": "ALL", "quoted": "UNCOVERED IN CORRECTIONS SECTION",
             "n_distinct_uncovered": len(uncovered), "n_in_exempt_corrections": len(only_corr)})
rows.append({"document": "ALL", "quoted": "UNCOVERED EVERYWHERE", "n_unsourced_mentions": per_doc,
             "n_distinct_uncovered": len(uncovered),
             "uncovered": ";".join(sorted(uncovered)[:40])})

print("\n2e. how well is each PASSING figure actually constrained? the average hides this")
if figrates:
    arr = np.array([r[0] for r in figrates])
    print(f"   {len(figrates)} passing figures; per-figure chance-match rate median {np.median(arr):.1%}, "
          f"mean {arr.mean():.1%}, worst {arr.max():.0%}")
    for thr, gloss in ((0.5, "a pass says almost nothing"), (0.25, "a pass is weak evidence"),
                       (0.10, "a pass is suggestive rather than conclusive")):
        print(f"   {int((arr > thr).sum()):4d} figures are above {thr:.0%} — {gloss}")
    print("   the five least constrained, which are the ones to quote more precisely or cite more narrowly:")
    for r in sorted(figrates, reverse=True)[:5]:
        print(f"      {r[0]:>4.0%}  {r[1].split('/')[-1]:28s} {r[2]:>9s} ({r[4]} dp) cites {r[3]}")
    rows.append({"document": "ALL", "quoted": "PER-FIGURE RATES", "n_checked": len(figrates),
                 "chance_match_cited_file": float(np.median(arr)),
                 "n_above_50pct": int((arr > 0.5).sum()), "n_above_10pct": int((arr > 0.10).sum())})
    print("   A figure above 50% is not verified in any useful sense: half the values at its precision are in")
    print("   its pool, so passing says almost nothing. Quoting it to more decimals is the cheapest fix, and")
    print("   citing one script rather than four is the next.")
    print("   ONE CLASS OF HIGH RATE IS NOT A DEFECT, and acting on it would be a mistake: a deliberately round")
    print("   DESIGN TARGET is not a measurement. docs/CHARGE_REVERSAL_DESIGN.md's −0.05 is a chosen number the")
    print("   experiment is sized against, and quoting it to four decimals would invent precision it does not")
    print("   have. The rate is a prompt to check whether a figure is a measurement quoted too coarsely, not an")
    print("   instruction to add digits.")

print("\n3. calibration: how often would a made-up coefficient match, i.e. how often would this check MISS a stale number?")
rng = np.random.default_rng(0)
RATE = {}
for dp in DPS:
    # rnd(), not round(), because a calibration that measures a DIFFERENT matching rule from the one the check
    # applies is not a calibration of that check. The first version of the five-decimal extension left round()
    # here after fixing the check itself, which would have reported the miss rate of a test nobody was running.
    per_file = [np.mean([rnd(v, dp) in rounded(pl, dp) for v in rng.uniform(-0.6, 0.6, 400)])
                for pl in file_pool.values()]
    cited = float(np.mean(per_file))
    glob_rate = float(np.mean([rnd(v, dp) in GLOBAL[dp] for v in rng.uniform(-0.6, 0.6, 1500)]))
    print(f"   {dp} dp: {cited:.1%} against a single cited file, {glob_rate:.0%} against every result file at once")
    RATE[dp] = cited
    rows.append({"document": "ALL", "quoted": f"CALIBRATION {dp}dp", "decimals": dp,
                 "chance_match_cited_file": cited, "chance_match_global": glob_rate})
print("   The second column is why the naive check is worthless and the first is why this one is not. Figures quoted to four")
print("   decimal places are the ones it genuinely constrains; at two it still misses about a fifth of stale numbers.")

# WHAT THE CALIBRATION MEANS FOR THIS CORPUS, not in the abstract. The per-precision chance-match rates above
# are easy to read as reassurance, and they are not: their weight depends entirely on how precisely the
# documents happen to quote. A figure given to two decimals is checked at a rate where a fifth of stale
# numbers would pass. So the useful statistic is the EXPECTED NUMBER of stale figures this check would pass
# given the precision actually used - a count, in the units of the thing being promised. It was computed
# after a stale pair of coefficients passed in one document and failed in another on the same day, purely
# because their pools differed: the clean run is evidence, not proof, and this says by how much.
from collections import Counter
prec = Counter(int(r["decimals"]) for r in rows
               if r.get("kind") != "sample size" and r.get("decimals") and r.get("found") is not None)
exp_missed = sum(n * RATE.get(dp, 0.0) for dp, n in prec.items())
tot_checked = sum(prec.values())
print("\n   what that means for THIS corpus, in figures rather than rates:")
for dp in sorted(prec):
    print(f"      {prec[dp]:4d} figures quoted to {dp} dp, chance match {RATE.get(dp, float('nan')):.1%}"
          f"  -> {prec[dp] * RATE.get(dp, 0.0):.1f} would pass if stale")
print(f"      so of {tot_checked} figures that pass, about {exp_missed:.0f} ({exp_missed / tot_checked:.0%}) "
      f"could be stale and pass undetected.")
print("   That is the number to quote when calling a run clean, and the way to lower it is to quote")
print("   coefficients to four decimals rather than two: the same check is nine times stronger there.")
rows.append({"document": "CALIBRATION", "quoted": "expected stale figures passed", "kind": "power",
             "n_checked": tot_checked, "n_untraceable": None,
             "chance_match_cited_file": exp_missed / tot_checked if tot_checked else None})

print(f"\n{'=' * 104}\n{len(all_missing)} coefficients across all documents cannot be found in the file their own claim cites:")
for doc, key, cites, ctx in all_missing:
    print(f"   {doc:34s} {key:>8s}  cites scripts/{cites}\n      {ctx}...")
if not all_missing:
    print("   (none)")
rows.append({"document": "ALL", "quoted": "SUMMARY", "n_untraceable": len(all_missing),
             "untraceable": ";".join(f"{d}:{k}" for d, k, _, _ in all_missing)})
pd.DataFrame(rows).to_csv(D / "document_consistency.csv", index=False)
print("\nwrote", D / "document_consistency.csv")
print("\nReading. A clean run means every coefficient asserted as current traces to a stored result. It does NOT mean the figure is")
print("the right one for the claim it is attached to, which no script judges, and it cannot see a non-numeric claim that is wrong.")
print("What it catches is the specific failure of correction 12: a number that outlived the analysis that produced it.")
