\name{sparse-correlations}
\alias{createSparseWithinOmicCorrelations}
\alias{createSparseCrossOmicCorrelations}
\alias{combineCorrelationNetworks}
\alias{sparseWithinOmicCorrelations}
\alias{sparseCrossOmicCorrelations}
\title{Sparse Prior-Guided Correlations}
\usage{
createSparseWithinOmicCorrelations(
  analysisData,
  knowledgeNetwork,
  assayNames = NULL,
  groupColumn = NULL,
  groupLevel = NULL,
  sampleIds = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  featureNameColumn = "featureName",
  resultName = NULL,
  storeResult = TRUE
)

createSparseCrossOmicCorrelations(
  analysisData,
  knowledgeNetwork,
  assayNames = NULL,
  groupColumn = NULL,
  groupLevel = NULL,
  sampleIds = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
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

sparseWithinOmicCorrelations(analysisData)

sparseCrossOmicCorrelations(analysisData)
}
\arguments{
  \item{analysisData}{A validated or validatable \code{MultiAssayExperiment}.}
  \item{knowledgeNetwork}{A standardized or standardizable prior network. Candidate edges must reference assay names present in \code{analysisData}.}
  \item{assayNames}{Optional character vector of assay names used to restrict candidate prior edges.}
  \item{groupColumn}{Optional sample metadata column used to select samples by group. Required when \code{sampleIds} is not supplied.}
  \item{groupLevel}{Single group label in \code{groupColumn}.}
  \item{sampleIds}{Optional explicit sample identifiers. Takes precedence over \code{groupColumn} and \code{groupLevel}.}
  \item{correlationMethod}{Correlation method. One of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.}
  \item{minimumAbsoluteCorrelation}{Numeric scalar or \code{NULL}. Edges with absolute correlation below this value are removed.}
  \item{adjustedPValueThreshold}{Numeric scalar or \code{NULL}. Edges with adjusted p-values above this value are removed.}
  \item{featureNameColumn}{Row metadata column containing display names for features.}
  \item{resultName}{Optional character scalar used when storing results.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return edge tables. If \code{TRUE}, store results in \code{metadata(analysisData)$cornetto} and return the updated \code{MultiAssayExperiment}.}
  \item{...}{One or more dense or sparse correlation edge tables, or lists of edge tables, to combine.}
  \item{removeDuplicates}{Logical scalar. If \code{TRUE}, duplicated interactions are removed using endpoint, edge type, source type, scope, method, group, and comparison metadata.}
}
\value{
\code{createSparseWithinOmicCorrelations()} and
\code{createSparseCrossOmicCorrelations()} return standardized edge
\code{DataFrame} objects when \code{storeResult = FALSE}; otherwise they
return an updated \code{MultiAssayExperiment}.

\code{combineCorrelationNetworks()} returns a standardized combined edge
\code{DataFrame} when \code{storeResult = FALSE}; otherwise it stores the
combined result in \code{metadata(analysisData)$cornetto$correlationNetworks}.

\code{sparseWithinOmicCorrelations()} and
\code{sparseCrossOmicCorrelations()} return stored sparse-correlation
results from the CorNetto metadata store.
}
\description{
Compute prior-restricted within-omic and cross-omic correlations and
combine dense and sparse dynamic network layers.
}
\details{
Sparse correlation functions do not test all possible feature pairs.
Instead, candidate pairs are taken from \code{knowledgeNetwork}. This is
useful when the number of possible pairs is large, when cross-omic all by
all testing is not biologically meaningful, or when the analysis should
be constrained to known or hypothesized interactions.

\code{createSparseWithinOmicCorrelations()} keeps only prior edges where
\code{fromAssayName == toAssayName}. \code{createSparseCrossOmicCorrelations()}
keeps only prior edges where \code{fromAssayName != toAssayName}. Both
functions compute pairwise correlations for the selected group, adjust
p-values with FDR correction, apply the requested correlation and
p-value thresholds, and return the same standardized edge schema used by
dense correlation functions.

If a prior network contains assays not present in \code{analysisData},
filter the prior network before calling sparse correlation functions.
Static priors from unmeasured assays can still be included later through
\code{createIntegratedNetwork()}.
}
\examples{
analysisData <- exampleAnalysisData()
knowledgeNetwork <- exampleKnowledgeNetwork()

measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))
measuredKnowledgeNetwork <- knowledgeNetwork[
    knowledgeNetwork$fromAssayName %in% measuredAssays &
        knowledgeNetwork$toAssayName %in% measuredAssays,
    ,
    drop = FALSE
]

sparseWithin <- createSparseWithinOmicCorrelations(
    analysisData = analysisData,
    knowledgeNetwork = measuredKnowledgeNetwork,
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)
head(sparseWithin)

sparseCross <- createSparseCrossOmicCorrelations(
    analysisData = analysisData,
    knowledgeNetwork = measuredKnowledgeNetwork,
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)
head(sparseCross)

dense <- createCorrelationNetwork(
    analysisData = analysisData,
    assayName = "protein",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)

combined <- combineCorrelationNetworks(dense, sparseWithin, sparseCross)
nrow(combined)

analysisData <- createSparseCrossOmicCorrelations(
    analysisData = analysisData,
    knowledgeNetwork = measuredKnowledgeNetwork,
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    resultName = "recoveredSparseCross",
    storeResult = TRUE
)
names(sparseCrossOmicCorrelations(analysisData))
}
