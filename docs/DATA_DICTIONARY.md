# Data dictionary (`data/templates/phenotype_template.csv`)

| Field | Meaning |
|---|---|
| strain_id | unique strain label |
| ncbi_assembly | exact matched GCA_/GCF_ accession; blank excludes the row from modelling |
| domain, phylum, class, family, genus, species | taxonomy; `family` is the blocking unit for CV, `class` for leave-one-lineage-out |
| ph_min, ph_opt, ph_max | pH limits and optimum, as reported |
| ph_min_bound, ph_max_bound | `measured` (growth failed beyond it), `lowest_tested`/`highest_tested` (censored), `extrapolated`, `unknown` |
| acid_identity | H2SO4, HCl, organic, mixed, buffered_unspecified, unspecified |
| organic_acid, organic_acid_mM | needed when acid_identity is organic or mixed |
| ph_measured_at_c | temperature of the pH reading (blank is flagged) |
| temp_opt_c, temp_min_c, temp_max_c | growth temperature |
| oxygen, medium, metabolism | culture context; metabolism e.g. Fe-oxidiser, S-oxidiser, heterotroph |
| metal_note | metal concentrations or resistance remarks |
| source, confidence | reference and high/medium/low |

Derived by `validate()`: `acid_class`, `acid_reach`, `ph_breadth`, `aH_min`, `ph_min_is_exact`.
Primary set = exact ph_min bound + ph_opt known + confidence high/medium + matched assembly.
