# CorNetto

CorNetto is a Bioconductor-style R package for knowledge-guided
multi-omic correlation network analysis. It is designed for normalized
transcriptomic, proteomic, and metabolomic abundance data stored in a
`MultiAssayExperiment`.

The first package version focuses on:

- group-specific correlation networks
- prior-guided sparse multi-omic correlations
- two-group differential correlation testing
- rewiring scores for perturbed nodes
- permutation validation for sparse edges and rewiring scores
- pathway-focused subnetworks
- graph construction and Cytoscape-ready export

## Installation

You can install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("bradleyalexward/CorNetto")
```

CorNetto depends on Bioconductor packages. If installation reports missing
Bioconductor dependencies, install them first:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "DGCA",
  "MultiAssayExperiment",
  "qvalue",
  "S4Vectors",
  "SummarizedExperiment"
))
```

## Citation

If you use or build upon CorNetto, please cite:

Ward, B. (2026). CorNetto: Knowledge-Guided Multi-Omic Correlation Network
Analysis. R package version 0.99.0.
https://github.com/bradleyalexward/CorNetto

You can also retrieve the package citation from R:

```r
citation("CorNetto")
```

A machine-readable citation is available in `CITATION.cff`.

## License

CorNetto is released under the Artistic License 2.0.
