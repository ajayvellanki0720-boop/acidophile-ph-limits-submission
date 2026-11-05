# Provenance

## Source snapshot

- Working project: `/Users/ajayvellanki/Downloads/Claude/acidophile-ph-limits`
- Upstream repository: `https://github.com/ajayvellanki0720-boop/acidophile-ph-limits`
- Requested and verified commit: `d03e70d88af04fc194b25c9005638407e23b687f`
- Snapshot date: 2026-11-05
- Submission repository history: new history created from the verified file tree, not a fork with the exploratory commit history.

The local tracked tree was exactly at the requested commit. The submission snapshot was exported from that commit, so ignored raw data and the original `.git` history were not copied.

## Local post-commit additions

The local working folder contained 13 untracked files created after the requested commit. They are included because they document additional checking of minimum-pH source statements, but they are not represented as inputs to the committed manuscript analyses:

- `scripts/107a_fetch_ph_statements.py` through `scripts/107f_recover_missing_sources.py`
- `results/broad2/ph_min_bacdive.csv`
- `results/broad2/ph_min_bacdive_observations.csv`
- `results/broad2/ph_min_fulltext.csv`
- `results/broad2/ph_min_recovered.csv`
- `results/broad2/ph_min_sources_missing.csv`
- `results/broad2/ph_min_subset_sensitivity.csv`
- `results/broad2/ph_min_verification.csv`

These additions are enumerated in `MANIFEST.tsv` and covered by `CHECKSUMS.sha256`.

## Generated submission files

The final Word manuscript is generated from `scripts/105_submission_docx.py` and then formatted to match the supplied PeerJ-style reference while retaining continuous native Word line numbering. The manuscript reports the new submission repository URL in its Data Summary.
