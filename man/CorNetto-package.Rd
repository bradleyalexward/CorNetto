\name{CorNetto-package}
\alias{CorNetto-package}
\title{CorNetto: Knowledge-Guided Multi-Omic Correlation Network Analysis}
\description{
CorNetto provides a Bioconductor-oriented workflow for constructing
knowledge-guided multi-omic correlation networks from normalized assay
data stored in a \code{MultiAssayExperiment}.
}
\details{
CorNetto is designed for studies where normalized transcriptomic,
proteomic, metabolomic, or related abundance assays are available for a
shared set of samples. The package reuses Bioconductor infrastructure,
including \code{MultiAssayExperiment}, \code{SummarizedExperiment}, and
\code{S4Vectors::DataFrame}, and returns standardized edge tables for
network analysis.

The main workflow is:
\enumerate{
  \item create and validate a \code{MultiAssayExperiment} with
  \code{createAnalysisData()} and \code{validateAnalysisData()};
  \item read or validate prior-knowledge networks with
  \code{readKnowledgeNetwork()} or \code{validateKnowledgeNetwork()};
  \item build dense within-omic networks with
  \code{createCorrelationNetwork()} or \code{createCorrelationNetworks()};
  \item add prior-guided sparse within-omic or cross-omic correlations
  with \code{createSparseWithinOmicCorrelations()} and
  \code{createSparseCrossOmicCorrelations()};
  \item compare two groups with \code{testDifferentialCorrelation()} and
  \code{createDifferentialCorrelationNetwork()};
  \item rank rewired nodes with \code{calculateRewiringScores()};
  \item combine static and dynamic layers with
  \code{createIntegratedNetwork()};
  \item extract pathway- or seed-focused subnetworks with
  \code{filterNetworkByNodes()} or \code{createFocusedNetwork()};
  \item create graph objects and export node and edge tables with
  \code{createNetworkGraph()}, \code{prepareCytoscapeTables()}, and
  \code{writeNetworkTables()}.
}
}
\seealso{
\code{\link{createAnalysisData}},
\code{\link{readKnowledgeNetwork}},
\code{\link{createCorrelationNetwork}},
\code{\link{createSparseCrossOmicCorrelations}},
\code{\link{testDifferentialCorrelation}},
\code{\link{createIntegratedNetwork}},
\code{\link{createFocusedNetwork}},
\code{\link{createNetworkGraph}}
}
\keyword{package}
