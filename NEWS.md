# CorNetto 0.99.0

- Initial Bioconductor-oriented package scaffold.
- Added support for normalized multi-omic inputs through
  `MultiAssayExperiment`.
- Implemented knowledge-network validation, dense and sparse correlation
  workflows, differential correlation testing, rewiring scores, pathway
  subnetworks, graph creation, and Cytoscape-ready export.
- Added synthetic example assets, tests, and an introductory vignette.
- Replaced separate sparse within-omic and cross-omic constructors with
  `createSparseMultiOmicCorrelations()`.
- Added permutation validation for sparse correlation edges and
  differential rewiring scores.
- Adjusted differential-correlation p-values across the full tested
  edge universe before applying significance filters.
