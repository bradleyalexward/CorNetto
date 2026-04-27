\name{network-integration}
\alias{createIntegratedNetwork}
\alias{filterNetworkByNodes}
\alias{createFocusedNetwork}
\alias{integratedNetworks}
\title{Integrated and Focused Networks}
\usage{
createIntegratedNetwork(
  analysisData = NULL,
  knowledgeNetwork = NULL,
  correlationNetwork = NULL,
  sparseMultiOmicCorrelation = NULL,
  differentialCorrelationNetwork = NULL,
  includeReverseEdges = TRUE,
  resultName = "integratedNetwork",
  storeResult = FALSE
)

filterNetworkByNodes(networkEdgeTable, nodes, mode = c("either", "both"))

createFocusedNetwork(
  networkEdgeTable,
  seedNodes,
  neighborhoodOrder = 1L,
  mode = c("all", "out", "in"),
  dropIsolatedNodes = FALSE
)

integratedNetworks(analysisData)
}
\arguments{
  \item{analysisData}{Optional \code{MultiAssayExperiment}. When supplied, stored CorNetto result layers can be retrieved automatically and integrated results can be stored back into metadata.}
  \item{knowledgeNetwork}{Optional standardized knowledge-network edge table or list of edge tables.}
  \item{correlationNetwork}{Optional dense correlation edge table or list of edge tables.}
  \item{sparseMultiOmicCorrelation}{Optional sparse multi-omic correlation edge table or list of edge tables.}
  \item{differentialCorrelationNetwork}{Optional differential-correlation edge table or list of edge tables.}
  \item{includeReverseEdges}{Logical scalar. If \code{TRUE}, undirected edges are duplicated in the reverse direction to support graph traversal and neighbor expansion.}
  \item{resultName}{Character scalar used when storing the integrated network.}
  \item{storeResult}{Logical scalar. If \code{FALSE}, return an edge table. If \code{TRUE}, return an updated \code{MultiAssayExperiment}.}
  \item{networkEdgeTable}{A standardized CorNetto edge table.}
  \item{nodes}{Character vector of node feature identifiers used by \code{filterNetworkByNodes()}.}
  \item{mode}{For \code{filterNetworkByNodes()}, either \code{"either"} to keep edges touching any selected node or \code{"both"} to keep edges where both endpoints are selected. For \code{createFocusedNetwork()}, neighborhood direction passed to \code{igraph::ego()}: \code{"all"}, \code{"out"}, or \code{"in"}.}
  \item{seedNodes}{Character vector of seed-node feature identifiers used to build a focused neighborhood.}
  \item{neighborhoodOrder}{Integer scalar. Number of graph steps to expand from \code{seedNodes}. Use \code{0} to keep only seed nodes.}
  \item{dropIsolatedNodes}{Logical scalar. If \code{TRUE}, nodes that are retained by the seed query but have no returned edges are removed from stored node metadata.}
}
\value{
\code{createIntegratedNetwork()}, \code{filterNetworkByNodes()}, and
\code{createFocusedNetwork()} return standardized edge \code{DataFrame}
objects. These edge tables may carry a stored node table in metadata so
isolated nodes can be preserved for graph creation and export.

When \code{createIntegratedNetwork(storeResult = TRUE)}, an updated
\code{MultiAssayExperiment} is returned.

\code{integratedNetworks()} returns stored integrated-network results
from \code{metadata(analysisData)$cornetto$integratedNetworks}.
}
\description{
Merge static prior knowledge and dynamic correlation layers, then filter
or expand the resulting network around node sets such as pathway seeds.
}
\details{
\code{createIntegratedNetwork()} row-binds any supplied network layers
after coercing them to the CorNetto standardized edge schema. Layers can
be supplied directly, retrieved from a \code{MultiAssayExperiment}, or
mixed. Static prior edges and dynamic correlation edges therefore remain
distinguishable through \code{sourceType}, \code{edgeType},
\code{correlationScope}, \code{correlationMethod}, \code{groupName}, and
\code{comparisonName}.

\code{filterNetworkByNodes()} performs direct endpoint filtering. Use
\code{mode = "either"} to keep first-contact edges around a node set, and
\code{mode = "both"} to keep only edges internal to the node set.

\code{createFocusedNetwork()} constructs a graph from the edge table,
finds all nodes within \code{neighborhoodOrder} steps of the seed nodes,
and returns the induced edge table. Node metadata is retained so a query
with \code{neighborhoodOrder = 0} can still return isolated seed nodes
for export.
}
\examples{
analysisData <- exampleAnalysisData()
knowledgeNetwork <- exampleKnowledgeNetwork()
seedNodes <- exampleSeedNodes()

dense <- createCorrelationNetwork(
    analysisData = analysisData,
    assayName = "protein",
    groupColumn = "clinicalGroup",
    groupLevel = "Recovered",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    storeResult = FALSE
)

integrated <- createIntegratedNetwork(
    knowledgeNetwork = knowledgeNetwork,
    correlationNetwork = dense,
    includeReverseEdges = TRUE,
    storeResult = FALSE
)
nrow(integrated)

directSubnetwork <- filterNetworkByNodes(
    networkEdgeTable = integrated,
    nodes = seedNodes,
    mode = "either"
)
nrow(directSubnetwork)

focusedNetwork <- createFocusedNetwork(
    networkEdgeTable = integrated,
    seedNodes = seedNodes,
    neighborhoodOrder = 1
)
nrow(focusedNetwork)

analysisData <- createIntegratedNetwork(
    analysisData = analysisData,
    knowledgeNetwork = knowledgeNetwork,
    correlationNetwork = dense,
    resultName = "exampleIntegrated",
    storeResult = TRUE
)
names(integratedNetworks(analysisData))
}
