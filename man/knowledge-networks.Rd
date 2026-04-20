\name{knowledge-networks}
\alias{readKnowledgeNetwork}
\alias{validateKnowledgeNetwork}
\alias{combineKnowledgeNetworks}
\alias{knowledgeNetworks}
\title{Prior-Knowledge Network Utilities}
\usage{
validateKnowledgeNetwork(
  knowledgeNetwork,
  analysisData = NULL,
  resultName = "knowledgeNetwork",
  storeResult = FALSE,
  edgeType = "association",
  edgeDirection = "undirected",
  knowledgeSource = "userSupplied"
)

readKnowledgeNetwork(
  filePath,
  columnMapping = NULL,
  delimiter = NULL,
  analysisData = NULL,
  resultName = "knowledgeNetwork",
  storeResult = FALSE,
  edgeType = "association",
  edgeDirection = "undirected",
  knowledgeSource = "userSupplied"
)

combineKnowledgeNetworks(
  ...,
  removeDuplicates = TRUE,
  analysisData = NULL,
  resultName = "combinedKnowledgeNetwork",
  storeResult = FALSE
)

knowledgeNetworks(analysisData)
}
\arguments{
  \item{knowledgeNetwork}{A data frame-like prior network. At minimum it must contain \code{fromFeatureIdentifier}, \code{toFeatureIdentifier}, \code{fromAssayName}, and \code{toAssayName}, unless those columns are created through \code{columnMapping} in \code{readKnowledgeNetwork()}.}
  \item{analysisData}{Optional \code{MultiAssayExperiment}. When \code{storeResult = TRUE}, the standardized network is stored in \code{metadata(analysisData)$cornetto}.}
  \item{resultName}{Character scalar used as the storage name when \code{storeResult = TRUE}.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return the standardized network. If \code{TRUE}, return an updated \code{MultiAssayExperiment}.}
  \item{edgeType}{Default edge type used when the input lacks an \code{edgeType} column.}
  \item{edgeDirection}{Default edge direction used when the input lacks an \code{edgeDirection} column. Typical values are \code{"undirected"}, \code{"directed"}, \code{"activation"}, or \code{"repression"}.}
  \item{knowledgeSource}{Default source label used when the input lacks a \code{knowledgeSource} column.}
  \item{filePath}{Character scalar. Path to a delimited prior-network file.}
  \item{columnMapping}{Optional named character vector. Names are standard CorNetto column names and values are column names in the input file.}
  \item{delimiter}{Optional delimiter override. When \code{NULL}, the delimiter is inferred from the file extension.}
  \item{...}{One or more prior-network objects, or lists of prior-network objects, to combine.}
  \item{removeDuplicates}{Logical scalar. If \code{TRUE}, duplicated rows are removed after combining networks.}
}
\value{
\code{validateKnowledgeNetwork()}, \code{readKnowledgeNetwork()}, and
\code{combineKnowledgeNetworks()} return a standardized edge
\code{DataFrame} when \code{storeResult = FALSE}. When
\code{storeResult = TRUE}, they return an updated
\code{MultiAssayExperiment}.

\code{knowledgeNetworks()} returns the list of stored knowledge-network
results from \code{metadata(analysisData)$cornetto$knowledgeNetworks}.
}
\description{
Read, validate, combine, and retrieve prior-knowledge networks used by
CorNetto to constrain sparse correlations and to build integrated
multi-omic networks.
}
\details{
The standardized CorNetto edge schema includes node identifiers and
display names for both endpoints, assay names for both endpoints,
edge-type metadata, source metadata, correlation fields, differential
correlation fields, edge weights, and directionality.

Only four columns are strictly required for a knowledge network:
\code{fromFeatureIdentifier}, \code{toFeatureIdentifier},
\code{fromAssayName}, and \code{toAssayName}. Optional columns such as
\code{fromFeatureName}, \code{toFeatureName}, \code{edgeType},
\code{edgeDirection}, \code{knowledgeSource}, \code{evidenceScore}, and
\code{isDirected} are preserved or filled from defaults.

Use \code{readKnowledgeNetwork()} when priors are stored on disk and need
column remapping. Use \code{validateKnowledgeNetwork()} when the object is
already in memory. Use \code{combineKnowledgeNetworks()} to merge
multiple prior sources, such as protein-protein, RNA-protein,
protein-metabolite, transcription-factor, or miRNA-target interactions.
}
\examples{
knowledgeFile <- system.file(
    "extdata",
    "example_knowledge_network.csv",
    package = "CorNetto",
    mustWork = TRUE
)

knowledgeNetwork <- readKnowledgeNetwork(knowledgeFile)
head(knowledgeNetwork)

inMemoryNetwork <- data.frame(
    fromFeatureIdentifier = "P1",
    toFeatureIdentifier = "P2",
    fromAssayName = "protein",
    toAssayName = "protein",
    edgeType = "proteinProteinInteraction",
    edgeDirection = "undirected",
    knowledgeSource = "example",
    stringsAsFactors = FALSE
)
validatedNetwork <- validateKnowledgeNetwork(inMemoryNetwork)

combinedNetwork <- combineKnowledgeNetworks(
    knowledgeNetwork,
    validatedNetwork,
    removeDuplicates = TRUE
)
nrow(combinedNetwork)

analysisData <- exampleAnalysisData()
analysisData <- validateKnowledgeNetwork(
    knowledgeNetwork = knowledgeNetwork,
    analysisData = analysisData,
    storeResult = TRUE
)
names(knowledgeNetworks(analysisData))
}
