\name{differential-correlation}
\alias{testDifferentialCorrelation}
\alias{createDifferentialCorrelationNetwork}
\alias{calculateRewiringScores}
\alias{differentialCorrelationResults}
\alias{rewiringResults}
\title{Differential Correlation and Rewiring}
\usage{
testDifferentialCorrelation(
  analysisData,
  groupColumn,
  groupLevels,
  assayName = NULL,
  candidateEdgeTable = NULL,
  featureSubset = NULL,
  correlationMethod = c("pearson", "spearman"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  pAdjustMethod = "qvalue",
  featureNameColumn = "featureName",
  resultName = NULL,
  storeResult = TRUE
)

createDifferentialCorrelationNetwork(
  differentialCorrelationTable,
  differenceAdjustedPValueThreshold = 0.05,
  minimumAbsoluteCorrelation = 0.3,
  edgeWeightMethod = c("scaledAbsoluteZScore", "absoluteZScore", "signedZScore")
)

calculateRewiringScores(
  differentialCorrelationNetwork,
  analysisData = NULL,
  resultName = "rewiringScores",
  storeResult = FALSE
)

differentialCorrelationResults(analysisData)

rewiringResults(analysisData)
}
\arguments{
  \item{analysisData}{A validated or validatable \code{MultiAssayExperiment}.}
  \item{groupColumn}{Sample metadata column defining the two groups to compare.}
  \item{groupLevels}{Character vector of length two. The first value is treated as group 1 and the second as group 2.}
  \item{assayName}{Assay name used when testing all within-omic feature pairs. Required when \code{candidateEdgeTable} is \code{NULL}.}
  \item{candidateEdgeTable}{Optional standardized candidate edge table. When supplied, only these candidate pairs are tested, and pairs can be within-omic or cross-omic.}
  \item{featureSubset}{Optional character vector of feature identifiers to retain before all-pairs within-omic testing. Ignored when \code{candidateEdgeTable} is supplied.}
  \item{correlationMethod}{Correlation method for differential testing. One of \code{"pearson"} or \code{"spearman"}.}
  \item{minimumAbsoluteCorrelation}{Numeric scalar. A candidate edge is retained only when at least one group has absolute correlation greater than or equal to this value.}
  \item{adjustedPValueThreshold}{Numeric scalar. A candidate edge is retained only when at least one group has adjusted within-group p-value less than or equal to this value.}
  \item{pAdjustMethod}{Multiple-testing correction method. Use \code{"qvalue"} to match the thesis workflow when appropriate, or any method accepted by \code{stats::p.adjust()}.}
  \item{featureNameColumn}{Row metadata column containing display names for features.}
  \item{resultName}{Optional character scalar used when storing results.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return edge tables. If \code{TRUE}, store results in \code{metadata(analysisData)$cornetto} and return the updated \code{MultiAssayExperiment}.}
  \item{differentialCorrelationTable}{Output of \code{testDifferentialCorrelation()}.}
  \item{differenceAdjustedPValueThreshold}{Numeric scalar. Maximum adjusted differential-correlation p-value retained in the network.}
  \item{edgeWeightMethod}{How edge weights are derived from the differential correlation z-score.}
  \item{differentialCorrelationNetwork}{A standardized differential-correlation edge table, usually from \code{createDifferentialCorrelationNetwork()}.}
}
\value{
\code{testDifferentialCorrelation()} returns a standardized edge
\code{DataFrame} containing group-specific correlations, group-specific
p-values, adjusted p-values, z-score differences, and differential
p-values when \code{storeResult = FALSE}; otherwise it returns an updated
\code{MultiAssayExperiment}.

\code{createDifferentialCorrelationNetwork()} returns a standardized
edge \code{DataFrame} with rewiring-oriented edge labels and edge
weights.

\code{calculateRewiringScores()} returns a \code{DataFrame} with one row
per node and columns including \code{nodeKey}, \code{nodeIdentifier},
\code{nodeName}, \code{assayName}, \code{totalConnections},
\code{rawRewiringScore}, \code{rootMeanSquareRewiringScore},
\code{degreeBin}, \code{degreeMatchedZScore}, and
\code{degreeMatchedScaledScore}. When \code{storeResult = TRUE}, it
returns an updated \code{MultiAssayExperiment}.

\code{differentialCorrelationResults()} and \code{rewiringResults()}
return stored results from the CorNetto metadata store.
}
\description{
Test whether correlations differ between two sample groups, convert
significant differences into a weighted differential-correlation network,
and summarize node-level rewiring.
}
\details{
\code{testDifferentialCorrelation()} supports two modes. In all-pairs
mode, it tests every feature pair in \code{assayName}, optionally after
restricting features with \code{featureSubset}. In candidate-edge mode,
it tests only the pairs supplied in \code{candidateEdgeTable}; this mode
is useful for prior-guided within-omic or cross-omic differential
correlation testing.

The output is filtered in two stages. First, an edge must have sufficient
absolute correlation in at least one group. Second, an edge must have a
sufficiently small adjusted within-group p-value in at least one group.
The differential p-value is also adjusted and stored in
\code{adjustedPValue}.

\code{createDifferentialCorrelationNetwork()} filters on the adjusted
differential-correlation p-value and classifies each edge as
\code{gainInGroup1}, \code{lossInGroup1}, \code{signFlip},
\code{strengthenedPositive}, \code{weakenedPositive},
\code{strengthenedNegative}, \code{weakenedNegative},
\code{positiveStable}, or \code{negativeStable}.

\code{calculateRewiringScores()} computes node-level summaries from the
weighted differential network. The raw rewiring score is the L2 norm of
incident edge weights. The root-mean-square score normalizes by node
degree. The degree-matched score standardizes raw scores within degree
bins to reduce degree-driven ranking effects.
}
\examples{
analysisData <- exampleAnalysisData()

differentialResults <- testDifferentialCorrelation(
    analysisData = analysisData,
    groupColumn = "clinicalGroup",
    groupLevels = c("PASC", "Recovered"),
    assayName = "protein",
    correlationMethod = "pearson",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    pAdjustMethod = "fdr",
    storeResult = FALSE
)
head(differentialResults)

differentialNetwork <- createDifferentialCorrelationNetwork(
    differentialCorrelationTable = differentialResults,
    differenceAdjustedPValueThreshold = 1,
    minimumAbsoluteCorrelation = 0
)
head(differentialNetwork)

rewiringTable <- calculateRewiringScores(
    differentialCorrelationNetwork = differentialNetwork,
    storeResult = FALSE
)
head(rewiringTable)

analysisData <- testDifferentialCorrelation(
    analysisData = analysisData,
    groupColumn = "clinicalGroup",
    groupLevels = c("PASC", "Recovered"),
    assayName = "protein",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    pAdjustMethod = "fdr",
    resultName = "pascVsRecoveredProtein",
    storeResult = TRUE
)
names(differentialCorrelationResults(analysisData))
}
