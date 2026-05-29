# CorNetto

CorNetto is a Bioconductor-style R package for knowledge-guided
multi-omic correlation network analysis. It is designed for normalized
transcriptomic, proteomic, and metabolomic abundance data stored in a
`MultiAssayExperiment`.

The first package version focuses on:

- group-specific correlation networks
- prior-guided sparse multi-omic correlations
- pre-analysis summaries of sample counts, missingness, and feature counts
- two-group differential correlation testing
- rewiring scores for perturbed nodes
- permutation validation for rewiring scores
- pathway-focused subnetworks
- graph construction and Cytoscape-ready export

## Installation

Before Bioconductor acceptance, install the development version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("bradleyalexward/CorNetto")
```

## Citation

If you use or build upon CorNetto, please cite:

Ward, B., Belkhir, L., Elens, L., HYGIEIA Consortium (2026). CorNetto: Knowledge-Guided Multi-Omic Correlation Network
Analysis. R package version 0.99.1.
https://github.com/bradleyalexward/CorNetto

You can also retrieve the package citation from R:

```r
citation("CorNetto")
```

A machine-readable citation is available in `CITATION.cff`.

## License

CorNetto is released under the Artistic License 2.0.
