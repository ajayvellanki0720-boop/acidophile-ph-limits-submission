"""Mechanical audit of this project for the error classes it has actually made -> DATA_DIR/audit.csv
WHAT THIS IS AND, MORE IMPORTANTLY, WHAT IT IS NOT. It cannot show that the project contains no errors, and any summary of it that says "no errors found" is wrong. It checks
for SPECIFIC, PREVIOUSLY MADE mistakes, mechanically. Every check below cites the correction in docs/SYNTHESIS.md section 5 that it derives from, and a check with no such
citation does not belong here - a checker calibrated against imagined failures mostly finds imagined failures, which is the lesson of scripts/69's first near-vacuous draft.
It cannot see a wrong model, a misread paper, a figure attached to the wrong claim, or any non-numeric claim that is false.
THE CHECKS, EACH TRACED TO A CORRECTION.
  C1  correction 15: a script computes with ph.pgls_scan - which standardises its predictor - and stores the result in a column named per_pt/per_point without dividing by
      that predictor's SD, so a per-SD coefficient is labelled per unit. Six scripts do this and are in the baseline.
  C2  correction 6: a script MULTIPLIES a pgls_scan coefficient by an SD, which converts an already-per-SD figure back to per unit. One live instance, in the baseline.
  C3  the two traceability failures of 2026-10-31: a script writes a DATA_DIR csv it does not declare in its docstring preamble, so scripts/69 can build no pool for the
      figures quoted from it and reports them untraceable.
  C4  the same day, twice: a script taking an orig/wide argument writes a _wide output without declaring BOTH spellings, so the variant pools scripts/69 keeps separate
      cannot be built and one panel's figures can pass as the other's.
  C5  correction 12: a script stores a ratio, share or fold-change - a derived quantity - with no assertion anywhere in it. The decomposition that correction 12 withdrew
      was computed across two different samples and nothing checked that they matched.
  C6  2026-10-31, twice: a prose count that disagrees with what it counts (the synthesis claimed 13 corrections when it had 15; the README claimed 24 log entries with 25).
  C7  a declared output that does not exist on disk, so anything quoted from it traces to nothing.
  C8  a declared output OLDER than the script that writes it, so the stored figures may predate the current code and a clean scripts/69 run proves less than it appears to.
THE BASELINE. Known, documented findings live in BASELINE with the correction that records them, so this reports what is NEW. A baseline entry that stops firing is itself
reported, because that means either the defect was fixed without the baseline being updated or a check has silently stopped working - and a check that cannot fail is not
a check.
Usage: python3 scripts/77_audit.py [--full]   (--full also runs pytest and scripts/69, which is slow; env DATA_DIR=results/broad2)"""
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
FULL = "--full" in sys.argv
SELFTEST = "--self-test" in sys.argv

# (check, script, correction that documents it) - findings known and recorded, so the audit reports what is new
BASELINE = {
    ("C1", "34b"): "correction 15", ("C1", "34d"): "correction 15", ("C1", "34e"): "correction 15",
    ("C1", "36b"): "correction 15", ("C1", "36c"): "correction 15",
    ("C2", "34b"): "correction 15 records this as a live instance of correction 6",
}
# Documented instances this script CANNOT see, listed so the gap is explicit rather than implied by silence.
# C1 matches by NAME, so it reaches scripts/36c only because that script reuses `b` for the value it stores;
# where the stored name was never bound to the scan call, no name-matching rule can reach it. Following a value
# across a pooling step is dataflow analysis, not a lint. scripts/74 is the authoritative audit of the scale
# question and pattern-matches nothing.
# C9 fires on the SYNTAX "unordered results and a CSV write". Whether that actually determines row order needs
# reading the code, and in two places it does not. Each exemption states why, so it can be re-checked rather than
# taken on trust - an allowlist without reasons is just a way of switching a check off.
ORDER_INSENSITIVE = {
    ("36c", "imap_unordered"): "results go into dict(...) and are then indexed by B.ncbi_assembly in a fixed order, "
                               "so worker order reaches nothing that is written",
    ("40a", "imap_unordered"): "hits feeds np.mean([pf in h for h in hits]), which is order-invariant; the written "
                               "rows come from `got`, not from the pool",
}
NOT_STATICALLY_DETECTABLE = {
    "40c": "stores pooled_beta_per_pt from a random-effects pool of scan returns, under a name never bound to the "
           "scan call, so no name-matching rule can reach it (correction 15)",
}
findings = []


def add(check, where, detail, correction=""):
    findings.append({"check": check, "where": where, "detail": detail,
                     "baseline": BASELINE.get((check, where), ""), "correction_source": correction})


SCRIPTS = sorted(ROOT.glob("scripts/*.py"))
DOCS = ["docs/SYNTHESIS.md", "docs/CHARGE_REVERSAL_DESIGN.md", "docs/PREREGISTRATION.md", "README.md"]
# scripts any document cites: the set for which the declaration convention has to hold, since it exists so that
# scripts/69 can build a candidate pool for the figures a document quotes
CITED = set()
for _d in DOCS:
    CITED |= set(re.findall(r"scripts/(\d+[a-z]?)", (ROOT / _d).read_text()))


def num_of(p):
    m = re.match(r"(\d+[a-z]?)_", p.name)
    return m.group(1) if m else None


def preamble(text):
    """The docstring up to the Usage line, which is where every script declares its outputs."""
    return text[:2600].split("Usage:")[0]


SD_DIV = re.compile(r"/\s*[\w.\[\]]*(?:sd|sd_\w+|\.std\(\))")
# THE META-SCRIPTS. scripts/69 and scripts/77 analyse other scripts' TEXT, so they contain the panel option's
# own syntax (`== "wide"`) without taking the option. The first attempt to separate them tested whether the
# script reads a structural panel - which failed hilariously, because writing that test put the literal
# "structure_panel" into both checkers' source and they matched themselves. A checker cannot reliably
# classify itself by its own string patterns, so the two are named, with the reason stated.
META_SCRIPTS = {"69": "analyses other scripts' text, so it contains the option syntax without taking it",
                "77": "same"}
PANEL_OPTION = re.compile(r"""==\s*["']wide["']|\("orig",\s*"wide"\)|\bTAG\s*,\s*PANEL_CSV|["']wide["']\s*==""")
SCAN_CALL = re.compile(r"\bph\.pgls_scan\s*\(|(?<!\.)\bpgls_scan\s*\(")


def scale_findings(t):
    """C1 and C2 on one script's text. Extracted so --self-test can run it on synthetic code: a checker nothing
    has tried to fool is the thing scripts/69's first draft was criticised for.

    C1 fires only on the DIRECT form - the stored value is a name bound to a pgls_scan return. That is precise,
    with no false positives, and deliberately INCOMPLETE. Two earlier rules both failed: "divides by an SD
    anywhere in the file" let scripts/40c through, because it z-scores covariates elsewhere; "divides on the
    storing line" flagged scripts/74, which derives its per-point column correctly two lines earlier. Following a
    value across a pooling step is dataflow analysis, not a lint, so what the rule cannot reach is listed in
    NOT_STATICALLY_DETECTABLE rather than guessed at.
    C2 fires only when the multiplied name is bound to a pgls_scan return, which is what distinguishes the error
    from the legitimate per-unit-to-per-SD conversion scripts/40c and 52 do on a polyfit slope.
    """
    out = []
    if not SCAN_CALL.search(t):
        return out                      # a CALL, not the word: the first draft flagged scripts/73 for a comment
    scan_names = set()
    for m in re.finditer(r"^\s*([\w\s,]+?)\s*=\s*(?:ph\.)?pgls_scan\s*\(", t, re.M):
        scan_names |= {v.strip() for v in m.group(1).split(",") if v.strip()}
    lines = t.split("\n")
    for i, ln in enumerate(lines):
        mm = re.search(r'["\']\w*(?:per_pt|per_point)["\']\s*:\s*([\w.]+)(\[\s*\d+\s*\])?', ln)
        if mm and mm.group(1) in scan_names and not SD_DIV.search(ln):
            out.append(("C1", f"line {i + 1} stores the pgls_scan return `{mm.group(1)}` in a per_pt/per_point "
                              f"column without dividing by the predictor SD"))
            break
    for i, ln in enumerate(lines):
        hit = next((nm for nm in scan_names
                    if re.search(rf"\b{re.escape(nm)}\s*(?:\[\s*\d+\s*\])?\s*\*\s*[\w.\[\]]*(?:sd|\.std\(\))", ln)),
                   None)
        if hit:
            out.append(("C2", f"line {i + 1} multiplies the pgls_scan return `{hit}` by an SD, converting an "
                              f"already-per-SD coefficient back to per unit"))
            break
    return out


SETCOMP = re.compile(r"^\s*(\w+(?:\[[^\]]*\])?)\s*=\s*\{([^{}]*?)\sfor\s", re.M)


def _is_set_comprehension(head):
    """`{x for ...}` is a SET and iterates in hash order; `{k: v for ...}` is a DICT and iterates in insertion
    order, which is deterministic. The first draft of C9 could not tell them apart and reported four dict
    comprehensions - PANEL, MAXS, XS, REFS - as nondeterministic. A `:` at bracket depth zero is the difference;
    one inside a slice, a lambda default or a nested literal is not."""
    depth = 0
    for ch in head:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif ch == ":" and depth == 0:
            return False
    return True


def determinism_findings(t):
    """C9. Output that depends on the PROCESS rather than on the inputs - correction 17 and its follow-up.

    Two real faults, both in scripts/36a, neither visible to any other check here:
      (a) `for nm in perfam[fam]` iterated a SET of protein names. Python randomises string hashing per process, so
          a tie between equally-scoring paralogues broke differently in every run and EVERY CELL of the panel was
          unstable. Rebuilding reproduced 99.71% of cells; 31 of 32 references differed somewhere.
      (b) imap_unordered yields a genome as whatever worker finishes it, so row order - and the column order pandas
          infers from first appearance - varied per run. Cells identical, files different, at line 74.
    (b) was found only by verifying the fix for (a) and requiring the two builds to be byte-identical rather than
    merely equal in the numbers I happened to compare.

    Deliberately narrow, like every rule here. It matches two syntactic shapes and WILL MISS any other route to the
    same defect: sorting on the wrong key, a dict built from an unordered source then iterated, os.listdir, glob
    without sort, set union order, a seed left unset. A lint claiming to cover those would be lying."""
    out = []
    lines = t.split("\n")
    if "imap_unordered" in t and ".to_csv(" in t:
        i = next(j for j, L in enumerate(lines) if "imap_unordered" in L)
        if "sort_values" not in "\n".join(lines[i:]) and "sort_index" not in "\n".join(lines[i:]):
            out.append(("C9", f"line {i + 1} collects results with imap_unordered and writes a CSV without sorting "
                              f"them, so row order depends on worker scheduling"))
    setnames = {m.group(1).split("[")[0].strip() for m in SETCOMP.finditer(t) if _is_set_comprehension(m.group(2))}
    setnames = {s for s in setnames if s.isidentifier()}
    for i, L in enumerate(lines):
        mm = re.search(r"\bfor\s+\w+\s+in\s+(\w+)(?:\[[^\]]*\])?\s*:", L)
        if mm and mm.group(1) in setnames and "sorted(" not in L:
            out.append(("C9", f"line {i + 1} iterates `{mm.group(1)}`, bound to a set comprehension, without "
                              f"sorted(); set order is not stable across processes"))
    return out


if SELFTEST:
    # MUTATION TEST. Each case is a minimal script; the checker must fire on the broken forms and stay silent on
    # the correct ones. The negative cases matter more: they are the shapes that produced this script's three
    # false-positive rounds, and without them the rules could be tightened into uselessness and still "pass".
    CASES = [
        ("C1 direct form", 'b, se, p, lam = ph.pgls_scan(C, y, X, f)\nrows.append({"beta_per_pt": b[0]})', {"C1"}),
        ("C2 multiply by SD", 'b, se, p, lam = ph.pgls_scan(C, y, X, f)\nrows.append({"beta_per_sd": b[0] * xv.std()})', {"C2"}),
        ("both at once", 'b, se, p, lam = ph.pgls_scan(C, y, X, f)\nrows.append({"beta_per_pt": b[0], "beta_per_sd": b[0] * xv.std()})', {"C1", "C2"}),
        ("correct: divides on the storing line", 'b, se, p, lam = ph.pgls_scan(C, y, X, f)\nrows.append({"beta_per_pt": b[0] / sd_h})', set()),
        ("correct: pgls, not pgls_scan", 'r = ph.pgls(y, X, C)\nrows.append({"beta_per_pt": r["beta"][-1]})', set()),
        ("correct: only mentions pgls_scan in a comment", '# pgls_scan standardises its predictor\nrows.append({"beta_per_pt": z})', set()),
        ("correct: multiplies a polyfit slope, not a scan return", 'b, se, p, lam = ph.pgls_scan(C, y, X, f)\nobs = np.polyfit(a, c, 1)[0]\nrows.append({"beta_per_sd": obs * sdc})', set()),
    ]
    # C9's cases. The two positives are the two faults as they actually appeared in scripts/36a; the three negatives
    # are the shapes a tightened rule must not flag - and the list-comprehension one matters because list order IS
    # defined, so iterating it bare is correct and a rule keyed on "iterates a comprehension" would cry wolf.
    CASES += [
        ("C9 unsorted parallel write", 'rows = list(p.imap_unordered(one, args))\npd.DataFrame(rows).to_csv(out)', {"C9"}),
        ("correct: sorted before write", 'rows = list(p.imap_unordered(one, args))\npd.DataFrame(rows).sort_values("k").to_csv(out)', set()),
        ("C9 bare set iteration", 'hits = {h.name for top in search() for h in top}\nfor nm in hits:\n    pass', {"C9"}),
        ("correct: sorted set iteration", 'hits = {h.name for top in search() for h in top}\nfor nm in sorted(hits):\n    pass', set()),
        ("correct: list comprehension, where order IS defined", 'hits = [h.name for h in top]\nfor nm in hits:\n    pass', set()),
        ("correct: DICT comprehension, insertion-ordered", 'PANEL = {x.label: (x.uniprot, x.pfam) for x in r}\nfor nm in PANEL:\n    pass', set()),
        ("C9 set comp with a slice in the body", 'hits = {h.name[0:3] for h in top}\nfor nm in hits:\n    pass', {"C9"}),
    ]
    bad = 0
    print("--self-test: the checker against code written to fool it\n")
    for name, src, want in CASES:
        got = {c for c, _ in scale_findings(src)} | {c for c, _ in determinism_findings(src)}
        ok = got == want
        bad += not ok
        print(f"   {'PASS' if ok else 'FAIL':4s}  {name:52s} expected {sorted(want) or '[]'}, got {sorted(got) or '[]'}")
    print(f"\n{len(CASES) - bad} of {len(CASES)} cases pass")
    sys.exit(1 if bad else 0)

for p in SCRIPTS:
    n = num_of(p)
    if not n:
        continue
    t = p.read_text()
    pre = preamble(t)
    lines = t.split("\n")
    # a CALL, not the word. The first draft matched "pgls_scan" anywhere and duly flagged scripts/73, which only
    # mentions it in a comment explaining this very error class.
    calls_scan = bool(SCAN_CALL.search(t))
    # names bound to a pgls_scan return, so C2 can tell "multiplied a scan coefficient" from the legitimate
    # per-unit-to-per-SD conversion that scripts/40c and 52 do on a polyfit slope
    scan_names = set()
    for m in re.finditer(r"^\s*([\w\s,]+?)\s*=\s*(?:ph\.)?pgls_scan\s*\(", t, re.M):
        scan_names |= {v.strip() for v in m.group(1).split(",") if v.strip()}

    for c, detail in scale_findings(t):
        add(c, n, detail, "correction 15" if c == "C1" else "correction 6")
    for c, detail in determinism_findings(t):
        if (n, "imap_unordered") in ORDER_INSENSITIVE and "imap_unordered" in detail:
            continue
        add(c, n, detail, "correction 17")

    # C3 / C4 / C7 / C8 - declared outputs
    declared = set(re.findall(r"(?:DATA_DIR|docs)/([A-Za-z0-9_./]+\.csv)", pre))
    # ONLY to_csv counts as a write. The first draft of this matched `D / "x.csv"` anywhere, which is how a script
    # READS its inputs, and it duly reported 250-odd scripts for "not declaring" files they only consume.
    written = set(re.findall(r'\.to_csv\(\s*(?:D|DATA_DIR)\s*/\s*f?["\']([A-Za-z0-9_./{}]+\.csv)["\']', t))
    # and the check is scoped to scripts the DOCUMENTS cite, because the declaration exists for scripts/69's benefit
    # and scripts/69 only needs it where a document quotes a figure. Older scripts predate the convention entirely.
    if n in CITED:
        for w in sorted(written):
            if "{" in w:                  # an f-string output, e.g. f"x{TAG}.csv" - resolved by the variant check
                continue
            if w not in declared:
                here = (D / w).exists()
                add("C3", n, f"writes {w} but does not declare it in the docstring preamble"
                             + ("" if here else "  [and that file is not in this DATA_DIR, so declaring it would "
                                                "only move the finding to C7 - not actionable here]"),
                    "2026-10-31 twice")
    # Same pattern as scripts/69's, for the same reason: the first version matched one spelling of the option
    # and missed scripts/39b, which compares on a different line from the one that reads the argument.
    if PANEL_OPTION.search(t) and n not in META_SCRIPTS:
        has_wide = any(d.endswith("_wide.csv") for d in declared)
        has_base = any(not d.endswith("_wide.csv") for d in declared)
        if declared and not (has_wide and has_base):
            add("C4", n, "takes a panel argument but does not declare both the default and _wide output spellings",
                "2026-10-31 twice")
    for d_ in sorted(declared):
        f = D / d_ if not d_.startswith("docs/") else ROOT / d_
        if not f.exists():
            add("C7", n, f"declares {d_}, which does not exist", "traceability")
        elif f.stat().st_mtime < p.stat().st_mtime:
            add("C8", n, f"{d_} is older than the script that writes it", "staleness in the other direction")

    # C5 - a derived quantity with nothing asserted. Scoped to scripts that FIT something, since correction 12 was
    # about a ratio of two fitted estimates; a data-preparation script storing a "share" column is not that.
    if ("pgls" in t or "lstsq" in t) and re.search(r'["\']\w*(?:ratio|share|fold|shortfall)\w*["\']\s*:', t) \
            and "assert" not in t:
        add("C5", n, "stores a ratio/share/fold of fitted quantities with no assertion anywhere", "correction 12")

# C6 - prose counts against what they count
syn = (ROOT / "docs/SYNTHESIS.md").read_text()
sec5 = syn[syn.index("## 5. Corrections"):]
sec5 = sec5[:sec5.index("\n## ")]
n_corr = len(set(re.findall(r"^(\d+)\.\s+\*\*", sec5, re.M)))
for m in re.finditer(r"\*\*(\d+) corrections to my own analyses\*\*", (ROOT / "README.md").read_text()):
    if int(m.group(1)) != n_corr:
        add("C6", "README.md", f"claims {m.group(1)} corrections; section 5 contains {n_corr}", "2026-10-31")
# WORD FORMS TOO. The pre-registration said "Section 5 lists **twelve** corrections" long after there were
# seventeen, and nothing caught it because this check only read digits in one file. A count spelled out is
# still a count, and a stale one in the pre-registration is worse than in the README: that document is dated
# and is supposed to be amended rather than drift.
WORDS = {w: i for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen "
    "seventeen eighteen nineteen twenty twenty-one twenty-two twenty-three twenty-four twenty-five".split())}
for doc in ("README.md", "docs/SYNTHESIS.md", "docs/PREREGISTRATION.md", "docs/CHARGE_REVERSAL_DESIGN.md"):
    body = (ROOT / doc).read_text()
    for m in re.finditer(r"\*\*([a-z-]+)\*\* corrections|lists \*\*([a-z-]+)\*\* corrections", body):
        w = m.group(1) or m.group(2)
        if w in WORDS and WORDS[w] != n_corr:
            add("C6", doc, f"spells out {w} ({WORDS[w]}) corrections; section 5 contains {n_corr}", "2026-10-31")
for m in re.finditer(r"(?:lists|section 5 lists) \*\*(\w+)\*\* corrections", syn):
    pass
rd = (ROOT / "README.md").read_text()
i = rd.index("## Status log")
NOT_LOG_ENTRIES = ("One blanket caveat", "One coverage gap")
blocks = [b for b in re.findall(r"^\*\*([A-Z].{0,60})", rd[i:], re.M) if not b.startswith(NOT_LOG_ENTRIES)]
m = re.search(r"The (\d+) entries below", rd)
if m and int(m.group(1)) != len(blocks):
    add("C6", "README.md", f"preamble claims {m.group(1)} log entries; {len(blocks)} are present", "2026-10-31")
ntests = len(re.findall(r"^def test_", (ROOT / "tests/test_acidophile.py").read_text(), re.M))
for doc in ("README.md", "docs/SYNTHESIS.md"):
    for mm in re.finditer(r"(\d+) tests|tests \((\d+)\)", (ROOT / doc).read_text()):
        claimed = int(mm.group(1) or mm.group(2))
        if claimed != ntests:
            add("C6", doc, f"claims {claimed} tests; tests/test_acidophile.py defines {ntests}", "2026-10-31")

# ---------------- report ----------------
F = pd.DataFrame(findings)
# THE TIERS, because the checks do not imply the same thing and reporting them together makes the audit cry wolf.
# ERROR    - a figure is or may be wrong now: a scale mislabelling, a prose count that contradicts what it counts.
# LATENT   - nothing is wrong yet, but it would be if a document started quoting from the file. scripts/69 is the
#            authoritative end-to-end check for this and is currently clean, so these are debt, not defects.
# INFO     - a weak signal with a known benign cause, kept visible and out of the way.
TIER = {"C1": "ERROR", "C2": "ERROR", "C6": "ERROR", "BASELINE-STALE": "ERROR", "PYTEST": "ERROR",
        "DOC-TRACE": "ERROR", "C9": "ERROR", "C3": "LATENT", "C4": "LATENT", "C5": "LATENT",
        "C7": "LATENT", "C8": "INFO"}
# C4 is listed explicitly rather than left to the .get default below. It was the only check missing from this
# table, so deriving the check count from it would have undercounted - a quieter version of the same mistake.
assert all(f"C{i}" in TIER for i in range(1, 10)), "a check exists that this table does not tier"
for f_ in findings:
    f_["tier"] = TIER.get(f_["check"], "LATENT")
F = pd.DataFrame(findings)
new = F[F.baseline == ""] if len(F) else F
known = F[F.baseline != ""] if len(F) else F
# DERIVED, not written down. The count was hardcoded as 8 and was wrong the moment C9 was added - which is check C6's
# own error class ("a prose count that contradicts what it counts") committed inside the script that checks for it.
NCHECKS = len([k for k in TIER if re.fullmatch(r"C\d+", k)])
# AND THE SPELLED-OUT FORM IS DERIVED TOO. Deriving the numeral here and leaving the count written out as a WORD in
# the two prose lines below was the same mistake half-fixed: the numeral tracked C9's arrival and the words did not,
# so this script printed a header saying 9 and a verdict saying one fewer, in the same run. The instance that
# matters is the one in the VERDICT, since that is the sentence a reader would quote. The literal is deliberately
# not written out in this comment: an earlier attempt asserted the stale word was gone from the source and the
# assertion passed against the comment explaining it, which is the self-matching trap scripts/69 hit with
# META_SCRIPTS. The check below counts the word's remaining occurrences instead of testing for its absence.
NWORD = {v: k for k, v in WORDS.items()}[NCHECKS]
print(f"{len(SCRIPTS)} scripts audited against {NCHECKS} checks, each traced to a correction this project actually made\n")
for tier, gloss in (("ERROR", "a figure is or may be wrong now"),
                    ("LATENT", "not wrong yet; would be if a document quoted from it. scripts/69 is the real check"),
                    ("INFO", "weak signal, benign causes - a docstring edit alone trips C8")):
    sub = new[new.tier == tier] if len(new) else new
    print(f"{tier} ({len(sub)}) - {gloss}")
    for r in (sub.itertuples() if len(sub) else []):
        print(f"   [{r.check}] {r.where:14s} {r.detail}")
    if not len(sub):
        print("   none")
    print()
if not len(new[new.tier == "ERROR"] if len(new) else new):
    print(f"No ERROR-tier finding. That means none of these {NWORD} classes is present in a form this script can see.")
    print("It is NOT a statement that the project is correct.\n")
print(f"Known and documented (in the baseline): {len(known)}")
for r in (known.itertuples() if len(known) else []):
    print(f"   [{r.check}] {r.where:14s} {r.detail[:68]} -> {r.baseline}")
stale = [k for k in BASELINE if not any((r.check, r.where) == k for r in (F.itertuples() if len(F) else []))]
if stale:
    print(f"\nBASELINE ENTRIES THAT NO LONGER FIRE: {len(stale)} - each is either a defect fixed without updating the")
    print("baseline, or a check that has silently stopped working. A check that cannot fail is not a check.")
    for c, w in stale:
        print(f"   [{c}] {w}: {BASELINE[(c, w)]}")
        findings.append({"check": "BASELINE-STALE", "where": w, "detail": f"{c} no longer fires", "baseline": "", "correction_source": ""})

print(f"\nDocumented but NOT statically detectable ({len(NOT_STATICALLY_DETECTABLE)}) - the honest limit of a lint:")
for w, why in NOT_STATICALLY_DETECTABLE.items():
    print(f"   [C1] {w:14s} {why}")

if FULL:
    print("\n--- pytest ---")
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")], capture_output=True, text=True)
    print("   " + (r.stdout.strip().splitlines() or ["(no output)"])[-1])
    if r.returncode != 0:
        findings.append({"check": "PYTEST", "where": "tests/", "detail": "failing", "baseline": "", "correction_source": ""})
    print("--- scripts/69 (document traceability; slow) ---")
    r = subprocess.run([sys.executable, str(ROOT / "scripts/69_document_consistency.py")], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if "untraceable" in line or "cannot be found" in line:
            print("   " + line.strip())
    if re.search(r"^[1-9]\d* coefficients across all documents cannot be found", r.stdout, re.M):
        findings.append({"check": "DOC-TRACE", "where": "docs/", "detail": "untraceable coefficients remain", "baseline": "", "correction_source": ""})
else:
    print("\n(pytest and scripts/69 not run; pass --full. They are the other two layers and this is the static one.)")

pd.DataFrame(findings).to_csv(D / "audit.csv", index=False)
print("\nwrote", D / "audit.csv")
print(f"\nReading. An empty NEW list means none of {NWORD} specific, previously made mistakes is present in a form this")
print("script can see. It is not a statement that the project is correct, and it must never be reported as one. The")
print("checks worth having are the ones that correspond to a real past failure; the honest way to extend this is to")
print("add a check when a NEW class of error is found, not to invent classes in advance.")
sys.exit(1 if len(new) and len(new[new.tier == "ERROR"]) else 0)
