\name{input-data}
\alias{readAssayData}
\alias{createAnalysisData}
\alias{validateAnalysisData}
\alias{filterFeatures}
\alias{mapFeatureIdentifiers}
\alias{exampleAnalysisData}
\alias{exampleKnowledgeNetwork}
\alias{exampleSeedNodes}
\alias{corNettoResults}
\title{Input Data, Validation, and Example Helpers}
\usage{
readAssayData(
  filePath,
  assayName = NULL,
  featureIdentifierColumn = 1L,
  featureNameColumn = NULL,
  delimiter = NULL
)

createAnalysisData(
  assayList,
  sampleData,
  featureIdentifierColumns = NULL,
  featureNameColumns = NULL
)

validateAnalysisData(analysisData)

filterFeatures(
  analysisData,
  assayNames = NULL,
  featureSubset = NULL,
  topVariableFeatures = NULL,
  minimumVariance = NULL
)

mapFeatureIdentifiers(
  analysisData,
  assayName = NULL,
  featureIdentifiers = NULL,
  featureNameColumn = "featureName"
)

exampleAnalysisData()

exampleKnowledgeNetwork()

exampleSeedNodes()

corNettoResults(analysisData)
}
\arguments{
  \item{filePath}{Character scalar. Path to a comma-separated or tab-separated assay table.}
  \item{assayName}{Optional character scalar. In \code{readAssayData()} it is used as a human-readable assay label in messages; in \code{mapFeatureIdentifiers()} it selects one assay when supplied.}
  \item{featureIdentifierColumn}{Column name or numeric column index containing unique feature identifiers.}
  \item{featureNameColumn}{Optional column name or numeric column index containing display names. Missing or empty display names fall back to the feature identifier.}
  \item{delimiter}{Optional delimiter override. When \code{NULL}, the delimiter is inferred from the file extension.}
  \item{assayList}{Named list of normalized assay objects. Each element can be a matrix, data frame, or \code{SummarizedExperiment}. Rows are features and columns are samples.}
  \item{sampleData}{Sample metadata. Row names must be sample identifiers, or a \code{sampleId} column must be present.}
  \item{featureIdentifierColumns}{Optional scalar or named list giving feature-identifier columns for data-frame assays in \code{assayList}.}
  \item{featureNameColumns}{Optional scalar or named list giving feature-name columns for data-frame assays in \code{assayList}.}
  \item{analysisData}{A \code{MultiAssayExperiment} to validate, filter, inspect, or retrieve stored CorNetto results from.}
  \item{assayNames}{Optional character vector of assay names to filter. When \code{NULL}, all assays are considered.}
  \item{featureSubset}{Optional character vector, or named list of character vectors, giving feature identifiers to retain.}
  \item{topVariableFeatures}{Optional integer scalar, or named list of integer scalars, giving the number of most variable features to retain per assay.}
  \item{minimumVariance}{Optional numeric scalar, or named list of numeric scalars, giving the minimum row variance required for retention.}
  \item{featureIdentifiers}{Optional character vector of feature identifiers to subset the mapping returned by \code{mapFeatureIdentifiers()}.}
}
\value{
\code{readAssayData()} returns a \code{SummarizedExperiment} containing
one assay named \code{abundance} and row metadata columns
\code{featureIdentifier} and \code{featureName}.

\code{createAnalysisData()}, \code{validateAnalysisData()}, and
\code{filterFeatures()} return a validated \code{MultiAssayExperiment}.

\code{mapFeatureIdentifiers()} returns a \code{DataFrame} containing
assay names, feature identifiers, and feature display names.

\code{exampleAnalysisData()} returns a synthetic
\code{MultiAssayExperiment}; \code{exampleKnowledgeNetwork()} returns a
standardized prior-network \code{DataFrame}; and
\code{exampleSeedNodes()} returns a character vector of seed-node
identifiers.

\code{corNettoResults()} returns the CorNetto result store from
\code{metadata(analysisData)$cornetto}.
}
\description{
Functions for reading normalized abundance tables, constructing the main
CorNetto \code{MultiAssayExperiment}, validating assay and sample
metadata, reducing feature sets before network inference, retrieving
feature labels, and loading small synthetic example objects.
}
\details{
\code{readAssayData()} expects a rectangular abundance table where one
column contains feature identifiers and the remaining non-metadata
columns are sample abundances. Values are coerced to numeric.

\code{createAnalysisData()} converts each assay to a
\code{SummarizedExperiment}, stores all assays in a
\code{MultiAssayExperiment}, initializes \code{metadata(analysisData)$cornetto},
and then validates the result.

\code{validateAnalysisData()} checks that sample identifiers are unique,
assay matrices are numeric, feature identifiers and sample identifiers
are unique, and all assay samples are represented in \code{colData}.

\code{filterFeatures()} reduces one or more assays by explicit feature
identifier, minimum variance, or top variable feature count. When several
filters are supplied, they are applied sequentially.

\code{mapFeatureIdentifiers()} returns assay-specific identifier to
display-name mappings using row metadata, which is useful when preparing
plots, Cytoscape tables, or manually checking prior-network mappings.
}
\examples{
proteinFile <- system.file(
    "extdata",
    "example_protein_abundance.csv",
    package = "CorNetto",
    mustWork = TRUE
)

proteinExperiment <- readAssayData(
    filePath = proteinFile,
    featureIdentifierColumn = "featureIdentifier",
    featureNameColumn = "featureName"
)
proteinExperiment

analysisData <- exampleAnalysisData()
analysisData <- validateAnalysisData(analysisData)

filteredData <- filterFeatures(
    analysisData = analysisData,
    assayNames = "protein",
    topVariableFeatures = 3
)
filteredData

featureMap <- mapFeatureIdentifiers(
    analysisData = analysisData,
    assayName = "protein"
)
head(featureMap)

knowledgeNetwork <- exampleKnowledgeNetwork()
seedNodes <- exampleSeedNodes()
head(knowledgeNetwork)
seedNodes

names(corNettoResults(analysisData))
}
