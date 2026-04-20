\name{network-graph}
\alias{createNetworkGraph}
\alias{plotRewiringScores}
\alias{prepareCytoscapeTables}
\alias{writeNetworkTables}
\title{Graph Creation, Rewiring Plots, and Network Export}
\usage{
createNetworkGraph(networkEdgeTable, nodeTable = NULL, directed = NULL)

plotRewiringScores(
  rewiringTable,
  scoreColumn = "degreeMatchedScaledScore",
  topN = 20L,
  label = NULL,
  main = "Top rewired nodes"
)

prepareCytoscapeTables(networkEdgeTable, rewiringTable = NULL)

writeNetworkTables(
  networkTables,
  directoryPath = ".",
  prefix = "cornettoNetwork",
  fileFormat = c("tsv", "csv")
)
}
\arguments{
  \item{networkEdgeTable}{A standardized CorNetto edge table.}
  \item{nodeTable}{Optional node table. When \code{NULL}, stored node metadata is used when present; otherwise a node table is derived from edge endpoints.}
  \item{directed}{Optional logical scalar controlling graph direction. When \code{NULL}, the graph is directed if any edge is marked as directed.}
  \item{rewiringTable}{Output of \code{calculateRewiringScores()}, or another node-level table containing the selected score column.}
  \item{scoreColumn}{Character scalar naming the rewiring-score column to plot.}
  \item{topN}{Integer scalar. Number of top-scoring nodes to display.}
  \item{label}{Optional label source for \code{plotRewiringScores()}. Use \code{"identifier"}, \code{"name"}, or \code{"key"}. When \code{NULL}, \code{nodeKey} is used if present, otherwise \code{nodeIdentifier}.}
  \item{main}{Character scalar plot title.}
  \item{networkTables}{A named list containing \code{nodes} and \code{edges}, usually from \code{prepareCytoscapeTables()}.}
  \item{directoryPath}{Output directory for written node and edge tables.}
  \item{prefix}{File-name prefix used by \code{writeNetworkTables()}.}
  \item{fileFormat}{Output file format. Either \code{"tsv"} or \code{"csv"}.}
}
\value{
\code{createNetworkGraph()} returns an \code{igraph} object.

\code{plotRewiringScores()} returns the plotted subset invisibly.

\code{prepareCytoscapeTables()} returns a named list with
\code{nodes} and \code{edges} elements.

\code{writeNetworkTables()} writes node and edge tables and invisibly
returns a named character vector containing the file paths.
}
\description{
Convert CorNetto edge tables to graph objects, visualize rewiring-score
rankings, prepare Cytoscape-style tables, and write node and edge tables
to disk.
}
\details{
\code{createNetworkGraph()} uses standardized edge endpoints to create an
\code{igraph} object. If the edge table contains stored node metadata,
that metadata is used so isolated nodes can be retained.

\code{plotRewiringScores()} orders a rewiring table by the requested
score column, keeps the top \code{topN} nodes, and displays a base R
barplot.

\code{prepareCytoscapeTables()} creates a node table and an edge table
with \code{source} and \code{target} columns suitable for Cytoscape
import. If \code{rewiringTable} is supplied, node-level rewiring columns
are merged into the node table.

\code{writeNetworkTables()} is the only exported CorNetto function that
writes network outputs to disk. It creates \code{directoryPath} if needed
and writes one node table and one edge table.
}
\examples{
analysisData <- exampleAnalysisData()

differentialResults <- testDifferentialCorrelation(
    analysisData = analysisData,
    groupColumn = "clinicalGroup",
    groupLevels = c("PASC", "Recovered"),
    assayName = "protein",
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = 1,
    pAdjustMethod = "fdr",
    storeResult = FALSE
)

differentialNetwork <- createDifferentialCorrelationNetwork(
    differentialCorrelationTable = differentialResults,
    differenceAdjustedPValueThreshold = 1,
    minimumAbsoluteCorrelation = 0
)

rewiringTable <- calculateRewiringScores(differentialNetwork)
head(plotRewiringScores(
    rewiringTable = rewiringTable,
    scoreColumn = "rawRewiringScore",
    topN = 3,
    label = "identifier"
))

graphObject <- createNetworkGraph(differentialNetwork)
graphObject

networkTables <- prepareCytoscapeTables(
    networkEdgeTable = differentialNetwork,
    rewiringTable = rewiringTable
)
names(networkTables)

writtenFiles <- writeNetworkTables(
    networkTables = networkTables,
    directoryPath = tempdir(),
    prefix = "cornettoExample",
    fileFormat = "tsv"
)
basename(writtenFiles)
}
