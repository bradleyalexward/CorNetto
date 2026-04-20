\name{correlation-networks}
\alias{createCorrelationNetwork}
\alias{createCorrelationNetworks}
\alias{correlationNetworks}
\title{Dense Within-Omic Correlation Networks}
\usage{
createCorrelationNetwork(
  analysisData,
  assayName,
  groupColumn = NULL,
  groupLevel = NULL,
  sampleIds = NULL,
  featureSubset = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  pAdjustMethod = "fdr",
  featureNameColumn = "featureName",
  resultName = NULL,
  storeResult = TRUE
)

createCorrelationNetworks(
  analysisData,
  assayNames = NULL,
  groupColumn = NULL,
  groupLevels = NULL,
  sampleIdGroups = NULL,
  featureSubset = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  pAdjustMethod = "fdr",
  featureNameColumn = "featureName",
  storeResult = TRUE
)

correlationNetworks(analysisData)
}
\arguments{
  \item{analysisData}{A validated or validatable \code{MultiAssayExperiment}.}
  \item{assayName}{Character scalar naming the assay to analyse with \code{createCorrelationNetwork()}.}
  \item{assayNames}{Optional character vector of assay names to analyse with \code{createCorrelationNetworks()}. When \code{NULL}, all assays are analysed.}
  \item{groupColumn}{Optional sample metadata column used to select samples by group. Required when \code{sampleIds} or \code{sampleIdGroups} are not supplied.}
  \item{groupLevel}{Single group label in \code{groupColumn}.}
  \item{groupLevels}{Optional character vector of group labels to analyse. When \code{NULL}, all unique values in \code{groupColumn} are used.}
  \item{sampleIds}{Optional explicit sample identifiers for one group-specific network. Takes precedence over \code{groupColumn} and \code{groupLevel}.}
  \item{sampleIdGroups}{Optional named list where each element is a character vector of sample identifiers. This allows custom group definitions without adding a group column to \code{colData}.}
  \item{featureSubset}{Optional character vector, or named list of character vectors, giving feature identifiers to retain before computing correlations.}
  \item{correlationMethod}{Correlation method. One of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.}
  \item{minimumAbsoluteCorrelation}{Numeric scalar or \code{NULL}. Edges with absolute correlation below this value are removed.}
  \item{adjustedPValueThreshold}{Numeric scalar or \code{NULL}. Edges with adjusted p-values above this value are removed.}
  \item{pAdjustMethod}{Multiple-testing correction method passed to \code{stats::p.adjust()}.}
  \item{featureNameColumn}{Row metadata column containing display names for features.}
  \item{resultName}{Optional character scalar used when storing a single network result.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return edge tables. If \code{TRUE}, store results in \code{metadata(analysisData)$cornetto} and return the updated \code{MultiAssayExperiment}.}
}
\value{
\code{createCorrelationNetwork()} returns a standardized edge
\code{DataFrame} when \code{storeResult = FALSE}; otherwise it returns
an updated \code{MultiAssayExperiment}.

\code{createCorrelationNetworks()} returns a named list of standardized
edge \code{DataFrame} objects when \code{storeResult = FALSE};
otherwise it returns an updated \code{MultiAssayExperiment}.

\code{correlationNetworks()} returns stored dense correlation-network
results from \code{metadata(analysisData)$cornetto$correlationNetworks}.
}
\description{
Construct dense group-specific within-omic correlation networks from
normalized assay matrices.
}
\details{
\code{createCorrelationNetwork()} selects samples from one assay, computes
all pairwise feature-feature correlations, adjusts p-values across the
tested feature pairs, and returns a standardized edge table. Each edge
contains feature identifiers, feature names, assay names, correlation
value, p-value, adjusted p-value, sign-based edge direction, correlation
method, group name, and edge weight.

\code{createCorrelationNetworks()} repeats the same operation across
multiple assays and groups. Use this wrapper for complete study-group
network construction, and use \code{createCorrelationNetwork()} when
tuning parameters for one assay or one group.

The correlation backend is selected by \code{correlationMethod}. Pearson
correlations use matrix calculations and t-test p-values. Spearman and
Kendall correlations use pairwise \code{stats::cor.test()} calls.
}
\examples{
analysisData <- exampleAnalysisData()

proteinNetwork <- createCorrelationNetwork(
    analysisData = analysisData,
    assayName = "protein",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    correlationMethod = "pearson",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    pAdjustMethod = "fdr",
    storeResult = FALSE
)
head(proteinNetwork)

networkList <- createCorrelationNetworks(
    analysisData = analysisData,
    assayNames = c("protein", "metabolite"),
    groupColumn = "clinicalGroup",
    groupLevels = c("Recovered", "PASC"),
    correlationMethod = "spearman",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)
names(networkList)

analysisData <- createCorrelationNetwork(
    analysisData = analysisData,
    assayName = "protein",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    resultName = "recoveredProtein",
    storeResult = TRUE
)
names(correlationNetworks(analysisData))
}
