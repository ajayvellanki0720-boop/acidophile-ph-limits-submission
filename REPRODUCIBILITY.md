# Reproducibility instructions

## Prerequisites

Use `environment.yml` for the closest reconstruction of the recorded environment. The original release used Python 3.9.6 on Apple silicon macOS with NumPy 2.0.2, pandas 2.3.3, SciPy 1.13.1, Biopython 1.85, pyHMMER 0.12.3, and Matplotlib 3.9.4.

## Verify deposited files

Run `bash reproduce.sh verify`. This verifies every deposited checksum, runs the test suite, checks manuscript/result consistency, and runs the static audit. The historical release document says this tier is the meaningful reader check because it does not require the excluded raw inputs.

## Rebuild figures

Run `bash reproduce.sh figures`. This executes `scripts/100_figures.py`, which reads stored result tables and rewrites the numbered figure files in `results/broad2/figures/` as both PNG and SVG.

## Regenerate the Word manuscript source

Run `bash reproduce.sh manuscript`. The generator asserts the headline values against their stored results before writing `docs/submission/Periplasmic_surface_charge_manuscript.docx`; `scripts/108_peerj_submission.py` then applies the supplied reference manuscript's margins, typography, paragraph rhythm, and black heading styles, embeds all five figure panels, and retains native continuous Word line numbering.

## Full reconstruction

Full reconstruction requires the raw inputs listed in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). Numbered script headers state their prerequisites and output paths. `docs/SYNTHESIS.md` section 9 narrates the main execution chain. This tier was not rerun during submission packaging because the excluded source collection is approximately 17 GB and some curation decisions are manual.
