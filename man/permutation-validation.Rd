\name{permutation-validation}
\alias{validateSparseCorrelationEdges}
\alias{permuteDifferentialRewiring}
\alias{validationResults}
\title{Permutation Validation}
\usage{
validateSparseCorrelationEdges(
  analysisData,
  knowledgeNetwork,
  assayNames = NULL,
  correlationScope = c("all", "withinOmic", "crossOmic"),
  groupColumn = NULL,
  groupLevel = NULL,
  sampleIds = NULL,
  correlationMethod = c("pearson", "spearman", "kendall"),
  nPermutations = 1000L,
  nullModel = c("permuteEndpointSamples", "permuteAssaySampleLabels"),
  alternative = c("two.sided", "greater", "less"),
  minimumAbsoluteCorrelation = NULL,
  empiricalPValueThreshold = NULL,
  pAdjustMethod = "fdr",
  featureNameColumn = "featureName",
  seed = NULL,
  resultName = NULL,
  storeResult = FALSE
)

permuteDifferentialRewiring(
  analysisData,
  groupColumn,
  groupLevels,
  assayName = NULL,
  candidateEdgeTable = NULL,
  featureSubset = NULL,
  correlationMethod = c("pearson", "spearman"),
  minimumAbsoluteCorrelation = 0.3,
  adjustedPValueThreshold = 0.05,
  pAdjustMethod = "fdr",
  featureNameColumn = "featureName",
  differenceAdjustedPValueThreshold = 0.05,
  edgeWeightMethod = c("scaledAbsoluteZScore", "absoluteZScore", "signedZScore"),
  nPermutations = 1000L,
  blockColumn = NULL,
  seed = NULL,
  keepPermutationScores = FALSE,
  resultName = NULL,
  storeResult = FALSE
)

validationResults(analysisData)
}
\arguments{
  \item{analysisData}{A validated or validatable \code{MultiAssayExperiment}.}
  \item{knowledgeNetwork}{A standardized or standardizable prior network.}
  \item{assayNames}{Optional assay names used to restrict sparse candidate edges.}
  \item{correlationScope}{Which sparse candidate edges to validate: \code{"all"}, \code{"withinOmic"}, or \code{"crossOmic"}.}
  \item{groupColumn}{Sample metadata column defining the group or groups.}
  \item{groupLevel}{Single group label used for sparse-edge validation.}
  \item{sampleIds}{Optional explicit sample identifiers for sparse-edge validation.}
  \item{correlationMethod}{Correlation method.}
  \item{nPermutations}{Positive integer number of permutations.}
  \item{nullModel}{Sparse-edge null model. \code{"permuteEndpointSamples"} permutes the second endpoint separately for each edge. \code{"permuteAssaySampleLabels"} permutes sample labels once per target assay per permutation.}
  \item{alternative}{Empirical tail used for sparse-edge p-values.}
  \item{minimumAbsoluteCorrelation}{Optional minimum observed absolute correlation retained.}
  \item{empiricalPValueThreshold}{Optional maximum adjusted empirical sparse-edge p-value retained.}
  \item{pAdjustMethod}{Multiple-testing correction method.}
  \item{featureNameColumn}{Row metadata column containing display names.}
  \item{seed}{Optional random seed.}
  \item{resultName}{Optional storage name.}
  \item{storeResult}{Logical scalar. If \code{TRUE}, store the validation result in \code{metadata(analysisData)$cornetto$validationResults}.}
  \item{groupLevels}{Character vector of length two giving the groups for differential rewiring.}
  \item{assayName}{Assay name used when testing all within-omic differential pairs.}
  \item{candidateEdgeTable}{Optional candidate edge table for prior-guided differential testing.}
  \item{featureSubset}{Optional feature subset for all-pairs differential testing.}
  \item{adjustedPValueThreshold}{Maximum adjusted within-group p-value retained by \code{testDifferentialCorrelation()}.}
  \item{differenceAdjustedPValueThreshold}{Maximum adjusted differential-correlation p-value retained in the rewiring network.}
  \item{edgeWeightMethod}{How edge weights are derived from differential-correlation z-scores.}
  \item{blockColumn}{Optional sample metadata column. When supplied, group labels are permuted within each block.}
  \item{keepPermutationScores}{Whether to return the node-by-permutation rewiring score table.}
}
\value{
\code{validateSparseCorrelationEdges()} returns a validation
\code{DataFrame} with observed correlations, observed adjusted p-values,
empirical p-values, adjusted empirical p-values, and null summaries.

\code{permuteDifferentialRewiring()} returns a named list containing the
observed differential-correlation table, observed differential network,
rewiring validation table, and optionally the permutation score table.

\code{validationResults()} returns stored validation outputs from the
CorNetto metadata store.
}
\description{
Permutation-based validation for sparse correlation edges and
node-level differential rewiring scores.
}
\details{
\code{validateSparseCorrelationEdges()} keeps the prior edge universe
fixed and breaks sample alignment before recomputing each correlation.
This tests whether observed prior-supported correlations are stronger
than expected under the selected sample-alignment null model.

\code{permuteDifferentialRewiring()} permutes group labels, reruns
\code{testDifferentialCorrelation()}, rebuilds the differential network
with \code{createDifferentialCorrelationNetwork()}, recalculates rewiring
scores, and compares observed node scores with the permuted null
distribution. Missing nodes in a permuted network are assigned a rewiring
score of zero for empirical p-value calculation.
}
\examples{
analysisData <- exampleAnalysisData()
knowledgeNetwork <- exampleKnowledgeNetwork()

sparseValidation <- validateSparseCorrelationEdges(
    analysisData = analysisData,
    knowledgeNetwork = knowledgeNetwork,
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    nPermutations = 10,
    seed = 1
)
head(sparseValidation)

rewiringValidation <- permuteDifferentialRewiring(
    analysisData = analysisData,
    groupColumn = "clinicalGroup",
    groupLevels = c("PASC", "Recovered"),
    assayName = "protein",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    differenceAdjustedPValueThreshold = 1,
    nPermutations = 5,
    seed = 1
)
head(rewiringValidation$rewiringTable)
}
