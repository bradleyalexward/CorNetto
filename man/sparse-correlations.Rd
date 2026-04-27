\name{sparse-correlations}
\alias{createSparseMultiOmicCorrelations}
\alias{combineCorrelationNetworks}
\alias{sparseMultiOmicCorrelations}
\title{Sparse Prior-Guided Multi-Omic Correlations}
\usage{
createSparseMultiOmicCorrelations(
  analysisData,
  knowledgeNetwork,
  assayNames = NULL,
  correlationScope = c("all", "withinOmic", "crossOmic"),
  groupColumn = NULL,
  groupLevel = NULL,
  sampleIds = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  pAdjustMethod = "fdr",
  featureNameColumn = "featureName",
  resultName = NULL,
  storeResult = TRUE
)

combineCorrelationNetworks(
  ...,
  removeDuplicates = TRUE,
  analysisData = NULL,
  resultName = "combinedCorrelationNetwork",
  storeResult = FALSE
)

sparseMultiOmicCorrelations(analysisData)
}
\arguments{
  \item{analysisData}{A validated or validatable \code{MultiAssayExperiment}.}
  \item{knowledgeNetwork}{A standardized or standardizable prior network. Candidate edges must reference assay names present in \code{analysisData}; unmeasured assay pairs are ignored by \code{createSparseMultiOmicCorrelations()}.}
  \item{assayNames}{Optional character vector of assay names used to restrict candidate prior edges.}
  \item{correlationScope}{Which candidate edges to test: \code{"all"}, \code{"withinOmic"}, or \code{"crossOmic"}.}
  \item{groupColumn}{Optional sample metadata column used to select samples by group. Required when \code{sampleIds} is not supplied.}
  \item{groupLevel}{Single group label in \code{groupColumn}.}
  \item{sampleIds}{Optional explicit sample identifiers. Takes precedence over \code{groupColumn} and \code{groupLevel}.}
  \item{correlationMethod}{Correlation method. One of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.}
  \item{minimumAbsoluteCorrelation}{Numeric scalar or \code{NULL}. Edges with absolute correlation below this value are removed.}
  \item{adjustedPValueThreshold}{Numeric scalar or \code{NULL}. Edges with adjusted p-values above this value are removed.}
  \item{pAdjustMethod}{Multiple-testing correction method. Use \code{"qvalue"} or any method accepted by \code{stats::p.adjust()}.}
  \item{featureNameColumn}{Row metadata column containing display names for features.}
  \item{resultName}{Optional character scalar used when storing results.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return edge tables. If \code{TRUE}, store results in \code{metadata(analysisData)$cornetto} and return the updated \code{MultiAssayExperiment}.}
  \item{...}{One or more dense or sparse correlation edge tables, or lists of edge tables, to combine.}
  \item{removeDuplicates}{Logical scalar. If \code{TRUE}, duplicated interactions are removed using endpoint, edge type, source type, scope, method, group, and comparison metadata.}
}
\value{
\code{createSparseMultiOmicCorrelations()} returns a standardized edge
\code{DataFrame} when \code{storeResult = FALSE}; otherwise it returns
an updated \code{MultiAssayExperiment}.

\code{combineCorrelationNetworks()} returns a standardized combined edge
\code{DataFrame} when \code{storeResult = FALSE}; otherwise it stores the
combined result in \code{metadata(analysisData)$cornetto$correlationNetworks}.

\code{sparseMultiOmicCorrelations()} returns stored sparse-correlation
results from the CorNetto metadata store.
}
\description{
Compute prior-restricted sparse correlations across within-omic and
cross-omic candidate edges, and combine dense and sparse dynamic network
layers.
}
\details{
Sparse correlation functions do not test all possible feature pairs.
Instead, candidate pairs are taken from \code{knowledgeNetwork}. This is
useful when the number of possible pairs is large, when cross-omic all by
all testing is not biologically meaningful, or when the analysis should
be constrained to known or hypothesized interactions.

\code{createSparseMultiOmicCorrelations()} uses
\code{correlationScope = "withinOmic"} to keep only prior edges where
\code{fromAssayName == toAssayName}, \code{correlationScope = "crossOmic"}
to keep only prior edges where \code{fromAssayName != toAssayName}, and
\code{correlationScope = "all"} to test both classes. The returned edge
table stores per-edge scopes as \code{"sparseWithinOmic"} or
\code{"sparseCrossOmic"}.

P-values are adjusted across the full set of successfully tested
candidate edges before correlation or adjusted-p-value filters are
applied.
}
\examples{
analysisData <- exampleAnalysisData()
knowledgeNetwork <- exampleKnowledgeNetwork()

sparseMulti <- createSparseMultiOmicCorrelations(
    analysisData = analysisData,
    knowledgeNetwork = knowledgeNetwork,
    correlationScope = "all",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)
head(sparseMulti)

dense <- createCorrelationNetwork(
    analysisData = analysisData,
    assayName = "protein",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)

combined <- combineCorrelationNetworks(dense, sparseMulti)
nrow(combined)

analysisData <- createSparseMultiOmicCorrelations(
    analysisData = analysisData,
    knowledgeNetwork = knowledgeNetwork,
    correlationScope = "crossOmic",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    resultName = "recoveredSparseMulti",
    storeResult = TRUE
)
names(sparseMultiOmicCorrelations(analysisData))
}
