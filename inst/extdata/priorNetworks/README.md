# Synthetic candidate networks

The CSV files in this directory contain arbitrary edges generated solely for
the package vignettes and tests. They are not extracted from an interaction
database and must not be interpreted as biological evidence.
They were created for CorNetto and are distributed under the package's
Artistic-2.0 license.

`synthetic_within_assay.csv` contains undirected within-assay edges.
`synthetic_cross_assay.csv` contains directed cross-assay edges. The direction
is illustrative, not causal. `synthetic_decoy_edges.csv` contains endpoints
that are intentionally absent from the packaged assays so the filtering step
in the COVID workflow has observable work to do.

Run `inst/scripts/make-synthetic-priors.R` from the package source tree to
rebuild all three files deterministically from the packaged COVID matrices.
