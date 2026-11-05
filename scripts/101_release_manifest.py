"""The release manifest: what version of what produced which result -> docs/RELEASE.md and DATA_DIR/release_manifest.csv
WHAT A REPRODUCIBILITY PACKAGE IS FOR, which decides what goes in it. It is not a gesture toward openness; it is the
answer to one question a reader will actually ask - "if I ran this, would I get the same numbers, and how would I
know?" So the manifest records three things and refuses to record a fourth. It records the ENVIRONMENT, because a
coefficient is a function of the library that computed it. It records a CHECKSUM per stored result, because the
documents cite those files by name and a reader needs to know the file they have is the file the claim was checked
against. It records the RUN ORDER, because 100 scripts with inter-dependencies are not reproducible by inspection.
WHAT IT REFUSES TO RECORD is a claim that the whole pipeline re-runs from nothing. It does not, on this machine, in
any reasonable time: the inputs include 12 GB of genomes and 4.3 GB of proteomes that are gitignored by project
policy, and the phenotype curation involved hand-reading primary literature that no script replays. A manifest that
implied otherwise would be the most consequential false claim in the repository, because it is the one a reader is
least likely to check. The manifest therefore states plainly which layers are reproducible from the repository alone
and which need re-fetched inputs.
THE THREE TIERS, derived from what is actually committed rather than asserted:
  TIER 1, from the repository alone: the documents, their traceability check, the static audit and the test suite.
          A reader can verify that every quoted coefficient matches a stored result, and that the library passes its
          mutation-checked guards, without downloading anything.
  TIER 2, needs the stored results (737 MB, committed): every analysis script that reads DATA_DIR and writes a table.
          These re-run in minutes to hours and reproduce the published numbers exactly.
  TIER 3, needs re-fetched raw inputs: panel construction, HMM detection, phylogeny. Documented, not rehearsed.
Usage: python3 scripts/101_release_manifest.py   (env DATA_DIR=results/broad2)"""
import hashlib
import os
import platform
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import pandas as pd

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
rows = []

print("1. environment")
import matplotlib
import numpy
import scipy
import Bio
try:
    import pyhmmer
    HM = pyhmmer.__version__
except ImportError:
    HM = "not installed"
ENV = [("python", sys.version.split()[0]), ("platform", platform.platform()),
       ("numpy", numpy.__version__), ("pandas", pd.__version__), ("scipy", scipy.__version__),
       ("biopython", Bio.__version__), ("pyhmmer", HM), ("matplotlib", matplotlib.__version__)]
for k, v in ENV:
    print(f"   {k:12s} {v}")
    rows.append({"part": "environment", "name": k, "value": v})
try:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip())
except Exception:
    commit, dirty = "unknown", True
print(f"   commit       {commit[:12]}{'  (WORKING TREE DIRTY - this manifest describes uncommitted state)' if dirty else ''}")
rows.append({"part": "environment", "name": "commit", "value": commit, "clean": not dirty})

print("\n2. stored results, with checksums")
# Only the files the documents actually cite are checksummed. A manifest of every file in results/ would be longer
# and less useful: what a reader needs is the provenance of the tables the claims point at.
cited = set()
for doc in ("docs/SYNTHESIS.md", "docs/MANUSCRIPT.md", "docs/PREREGISTRATION.md", "docs/CHARGE_REVERSAL_DESIGN.md",
            "docs/DATABASES.md", "docs/EXPERIMENTAL_EVIDENCE.md", "README.md"):
    txt = (ROOT / doc).read_text()
    for m in re.finditer(r"`scripts/(\d+[a-z]?)(?:\s+(orig|wide))?`", txt):
        cited.add((m.group(1), m.group(2)))
decl = {}
for f in sorted((ROOT / "scripts").glob("*.py")):
    num = f.name.split("_")[0]
    head = f.read_text().split('"""')[1] if '"""' in f.read_text() else ""
    for out in re.findall(r"DATA_DIR/([A-Za-z0-9_./]+\.csv)", head.split("Usage:")[0]):
        decl.setdefault(num, set()).add(out)
    for out in re.findall(r"\bresults/([A-Za-z0-9_./]+\.csv)", head.split("Usage:")[0]):
        decl.setdefault(num, set()).add(out.split("/")[-1])
want = sorted({o for num, _ in cited for o in decl.get(num, ())})
n_ok = n_missing = 0
for o in want:
    fp = D / o
    if not fp.exists():
        rows.append({"part": "result", "name": o, "value": "MISSING", "bytes": 0})
        n_missing += 1
        continue
    h = hashlib.sha256(fp.read_bytes()).hexdigest()[:16]
    rows.append({"part": "result", "name": o, "value": h, "bytes": fp.stat().st_size})
    n_ok += 1
print(f"   {len(cited)} distinct script citations across the documents; {len(want)} declared result files")
print(f"   {n_ok} present and checksummed, {n_missing} declared but absent from this DATA_DIR")

print("\n3. the three tiers, and what each needs")
TIERS = [("1 - repository alone", "documents, traceability check, static audit, test suite",
          "python3 scripts/69_document_consistency.py; python3 scripts/77_audit.py; python3 -m pytest -q", "nothing"),
         ("2 - stored results", "every analysis script reading DATA_DIR",
          "python3 scripts/<n>_<name>.py [orig|wide]", f"results/broad2 ({_sz(D) if (_sz := (lambda p: f'{sum(f.stat().st_size for f in p.rglob(chr(42)) if f.is_file()) / 1e6:.0f} MB')) else ''})"),
         ("3 - raw inputs", "panel construction, HMM detection, phylogeny",
          "see docs/RELEASE.md", "12 GB genomes + 4.3 GB proteomes, gitignored, re-fetch documented")]
for t, what, how, needs in TIERS:
    print(f"   TIER {t:22s} {what}")
    print(f"        needs: {needs}")
    rows.append({"part": "tier", "name": t, "value": what, "how": how, "needs": needs})

man = pd.DataFrame(rows)
man.to_csv(D / "release_manifest.csv", index=False)
print(f"\nwrote {D / 'release_manifest.csv'}")

REL = ROOT / "docs/RELEASE.md"
env_tbl = "\n".join(f"| `{k}` | {v} |" for k, v in ENV)
res_tbl = "\n".join(f"| `{r['name']}` | `{r['value']}` | {r['bytes'] / 1000:.0f} kB |"
                    for r in rows if r["part"] == "result" and r["value"] != "MISSING")
REL.write_text(f"""# Reproducibility package

Generated by `scripts/101_release_manifest.py` on {date.today().isoformat()} at commit `{commit[:12]}`.
{'' if not dirty else chr(10) + '> **The working tree was dirty when this was generated**, so the checksums below describe uncommitted state. Regenerate after committing before releasing.' + chr(10)}
**What is and is not reproducible, stated before anything else.** The full pipeline does **not** re-run from this
repository. Its inputs include 12 GB of genome assemblies and 4.3 GB of proteomes that are gitignored by policy, and
the phenotype curation involved hand-reading primary literature that no script replays. A package implying otherwise
would be the most consequential false claim here, because it is the one a reader is least likely to check. What *is*
reproducible is set out in tiers below, and the first tier needs nothing but this repository.

## Environment

These are the versions that produced the committed results. The estimator is sensitive to none of them in principle
and was reproduced independently under the same ones (`scripts/95`), which is a weaker guarantee than cross-version
testing and is stated as such.

| Component | Version |
|---|---|
{env_tbl}

## Tiers

| Tier | Covers | Needs | How |
|---|---|---|---|
| **1** | the documents, their traceability check, the static audit, the test suite | nothing beyond this repository | `python3 scripts/69_document_consistency.py`<br>`python3 scripts/77_audit.py`<br>`python3 -m pytest -q` |
| **2** | every analysis script that reads `DATA_DIR` and writes a table | the committed `results/broad2` | `python3 scripts/<n>_<name>.py [orig\\|wide]` |
| **3** | structural panel, HMM detection, phylogeny | re-fetched raw inputs (see below) | documented, not rehearsed |

**A tier-1 run is the meaningful check for a reader.** It verifies that every signed coefficient in every document
matches the stored result the same claim cites, that none of nine previously-made error classes is mechanically
present, and that the library passes 26 tests including mutation-checked guards on the per-SD/per-unit identity.
It takes about fifteen seconds.

## Re-fetching the tier-3 inputs

Each script that needs a gitignored input states its own fetch command in a header comment. The two added most
recently:

```
curl -O https://www.bork.embl.de/~robbani/metatraits/ncbi_species_summary_no_predictions.tsv.gz
gzcat ncbi_species_summary_no_predictions.tsv.gz | awk -F'\\t' 'NR==1 || $3=="pH minimum"' > data/raw/metatraits/ncbi_species_ph_minimum_no_predictions.tsv
gzcat ncbi_species_summary_no_predictions.tsv.gz | awk -F'\\t' 'NR==1 || $3 ~ /^salinity/' > data/raw/metatraits/ncbi_species_salinity_no_predictions.tsv
```

```
mkdir -p data/raw/sensor_hmm
for pf in PF08918 PF16527 PF00497 PF00672; do
  curl -sSL -o data/raw/sensor_hmm/$pf.hmm.gz "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/$pf/?annotation=hmm"
  gunzip -f data/raw/sensor_hmm/$pf.hmm.gz
done
```

Genome assemblies and proteomes come from NCBI by the accessions in `data/curated/accession_used.csv`.

## Checksums of the cited result files

SHA-256, first 16 hex characters. These are the tables the documents' claims point at; a reader checking a quoted
coefficient should confirm the file they have is the file it was checked against.

| File | SHA-256 (16) | Size |
|---|---|---|
{res_tbl}

## Order

Scripts are numbered in execution order and `docs/SYNTHESIS.md` §9 narrates the main chain. Environment variables
`DATA_DIR` and `PHYLO_DIR` select the dataset; `PANEL_CSV` selects the structural panel. Scripts taking an
`orig`/`wide` argument default to `wide`, which is the panel every current figure uses.

## Figures

`python3 scripts/100_figures.py` redraws all five panels from the stored results into `results/broad2/figures/`.
Each panel asserts its headline values against the manuscript's text before drawing, so a figure cannot drift from
the prose it illustrates.
""")
print(f"wrote {REL}")
