#' CorNetto: Knowledge-Guided Multi-Omic Correlation Network Analysis
#'
#' @description
#' \if{html}{\figure{logo.png}{options: style='float: right' alt='' width='120'}}
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
#' - differential correlation, over all pairs in an assay or over
#'   prior-supported candidate edges
#' - rewiring score calculation
#' - permutation-based assessment of rewiring scores
#' - pathway-focused network extraction
#' - graph creation and Cytoscape-ready export
#'
#' The vignette (`vignette("CorNetto")`) opens with a diagram of how these
#' fit together.
#'
#' @references
#' Ward B, Balligand J-L, Bamps L, Bommer G, Cani PD, De Greef J,
#' Dewulf JP, Gatto L, Haufroid V, Kabamba B, Pyr dit Ruys S,
#' Vertommen D, Yombi JC, Belkhir L, Elens L (2026). Longitudinal
#' multi-omic network rewiring at the complement-coagulation interface
#' in post-acute sequelae of COVID-19 (PASC). \emph{medRxiv}.
#' \doi{10.64898/2026.07.14.26358048}
#'
#' @return No return value. This package-level help page documents the
#'   CorNetto package and points users to the main workflow functions.
#' @keywords internal
"_PACKAGE"
