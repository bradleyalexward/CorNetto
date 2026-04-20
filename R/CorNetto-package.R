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
#' - group-specific correlation networks
#' - prior-restricted sparse within-omic correlations
#' - prior-restricted sparse cross-omic correlations
#' - Pearson-based differential correlation analysis
#' - rewiring score calculation
#' - pathway-focused network extraction
#' - graph creation and Cytoscape-ready export
#'
#' @keywords internal
"_PACKAGE"
