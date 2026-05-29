#' CorNetto: Knowledge-Guided Multi-Omic Correlation Network Analysis
#'
#' CorNetto provides a Bioconductor-oriented workflow for constructing
#' knowledge-guided correlation networks from normalized multi-omic
#' abundance data. The package is designed around
#' `MultiAssayExperiment::MultiAssayExperiment()` so that multiple assays
#' and their shared sample metadata can be analysed consistently.
#'
#' The main workflow includes:
#'
#' - validation of normalized assay inputs
#' - pre-analysis summaries of missingness, sample counts, and feature counts
#' - sample and feature filtering before network construction
#' - group-specific correlation networks
#' - dense and prior-guided differential correlation analysis
#' - rewiring score calculation
#' - permutation validation for rewiring scores
#' - pathway-focused network extraction
#' - graph creation and Cytoscape-ready export
#'
#' @return No return value. This package-level help page documents the
#'   CorNetto package and points users to the main workflow functions.
#' @keywords internal
"_PACKAGE"
