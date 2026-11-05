# Data availability and raw-input retrieval

## Included

This repository includes the curated phenotype and accession tables used by the analyses, stored analysis intermediates needed for verification, result tables, phylogenetic outputs that fit within GitHub limits, and final figures.

The accession used for every genome is recorded in `data/curated/accession_used.csv`. The principal manuscript result tables and their role are documented in `results/broad2/release_manifest.csv` and `docs/DATA_DICTIONARY.md`.

## Intentionally excluded

The working project contains approximately 17 GB under `data/`, primarily downloaded genome assemblies, proteomes, structures, database exports, and other raw third-party inputs. These are excluded because they are re-fetchable, oversized for ordinary Git hosting, or not appropriate to redistribute. Executables under `tools/` are also excluded; install the named upstream tools instead.

## Fetch instructions

- NCBI assemblies and proteomes: retrieve the accessions listed in `data/curated/accession_used.csv`; see `scripts/02a_fetch_genomes.py` for the project retrieval implementation.
- BacDive records: see `scripts/00_fetch_bacdive.py` and `scripts/35a_fetch_bacdive_archaea.py`. BacDive access may require credentials and is subject to its terms.
- Europe PMC metadata and available full text: see `scripts/00c_fetch_europepmc.py`, `scripts/00c2_fetch_europepmc_by_year.py`, and `scripts/00h_fetch_fulltext.py`.
- Protein structures: retrieve the accessioned models from the AlphaFold Protein Structure Database; panel identifiers are recorded in `data/curated/structure_panel_refs.csv`.
- Pfam HMMs: retrieve the families used by the relevant scripts from InterPro/EMBL-EBI.
- metaTraits: the exact shell commands for the measured-only pH and salinity exports are recorded in `docs/RELEASE.md`.

Third-party services and database exports can change. Record download dates and hashes when reconstructing the raw-input tier.

## Manual curation boundary

Phenotype curation included hand-reading primary literature. Scripts support acquisition, parsing, and audits, but no automated workflow reproduces every human judgment. The curated deposited tables are therefore the authoritative analysis inputs.
