# Reproducibility instructions

## Prerequisites

Use `environment.yml` for the closest reconstruction of the recorded environment. The original release used Python 3.9.6 on Apple silicon macOS with NumPy 2.0.2, pandas 2.3.3, SciPy 1.13.1, Biopython 1.85, pyHMMER 0.12.3, and Matplotlib 3.9.4.

## Verify deposited files

Run `bash reproduce.sh verify`. This verifies every deposited checksum, runs the test suite, checks manuscript/result consistency, and runs the static audit. The historical release document says this tier is the meaningful reader check because it does not require the excluded raw inputs.

## Rebuild figures

Run `bash reproduce.sh figures`. This executes `scripts/100_figures.py`, which reads stored result tables and rewrites the numbered figure files in `results/broad2/figures/` as both PNG and SVG.

## Regenerate the Word manuscript source

Run `bash reproduce.sh figures` before `bash reproduce.sh manuscript`. The base generator asserts headline values against stored results; `scripts/108_peerj_submission.py` applies reference typography and `scripts/peerj_revision.py` applies the current structured abstract, methods, terminology, limitations and disclosures. The final Word copy retains continuous line numbering and includes legends without embedded figures. `scripts/109_submission_assets.py` exports four separate journal figures (including combined Figure 4A/B) and creates the 30-MB-limit-checked `Supplemental_Data_S1.zip` for direct journal upload. PNGs are rendered at 300 dpi; SVGs preserve vector graphics.

The primary stored panel and its later deterministic rebuild are distinct. The manuscript explicitly records that the canonical headline uses the original panel and the measured deterministic-rebuild sensitivity is +0.1858 versus +0.1848 per SD. No full-sample maximum-likelihood tree analysis was run. Historical source documents retain earlier language and are not the current submission manuscript.

## Full reconstruction

Full reconstruction requires the raw inputs listed in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). Numbered script headers state their prerequisites and output paths. `docs/SYNTHESIS.md` section 9 narrates the main execution chain. This tier was not rerun during submission packaging because the excluded source collection is approximately 17 GB and some curation decisions are manual.
