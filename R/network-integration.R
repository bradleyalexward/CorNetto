.nodeKey <- function(assayName, featureIdentifier) {
    paste(assayName, featureIdentifier, sep = "::")
}

.duplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    edgeData <- as.data.frame(edgeTable, stringsAsFactors = FALSE, check.names = FALSE)
    edgeIsDirected <- as.logical(edgeData$isDirected)
    edgeIsDirected[is.na(edgeIsDirected)] <- FALSE
    undirectedIndex <- !edgeIsDirected
    if (!any(undirectedIndex)) {
        return(edgeTable)
    }

    duplicateTable <- edgeData[undirectedIndex, , drop = FALSE]
    oldFromIdentifier <- duplicateTable$fromFeatureIdentifier
    oldFromName <- duplicateTable$fromFeatureName
    oldFromAssay <- duplicateTable$fromAssayName

    duplicateTable$fromFeatureIdentifier <- duplicateTable$toFeatureIdentifier
    duplicateTable$fromFeatureName <- duplicateTable$toFeatureName
    duplicateTable$fromAssayName <- duplicateTable$toAssayName
    duplicateTable$toFeatureIdentifier <- oldFromIdentifier
    duplicateTable$toFeatureName <- oldFromName
    duplicateTable$toAssayName <- oldFromAssay

    .coerceStandardEdgeTable(rbind(edgeData, duplicateTable), nodeTable = nodeTable)
}

.createNodeTable <- function(edgeTable, rewiringTable = NULL, nodeTable = NULL) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    if (is.null(nodeTable)) {
        nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    } else {
        nodeTable <- .coerceStandardNodeTable(nodeTable)
    }

    if (!is.null(rewiringTable)) {
        rewiringTable <- as.data.frame(rewiringTable, stringsAsFactors = FALSE, check.names = FALSE)
        joinColumn <- if ("nodeKey" %in% names(rewiringTable)) "nodeKey" else "nodeIdentifier"
        rewiringTable <- rewiringTable[
            ,
            setdiff(names(rewiringTable), setdiff(.standardNodeColumns(), joinColumn)),
            drop = FALSE
        ]
        nodeTable <- merge(
            x = as.data.frame(nodeTable, stringsAsFactors = FALSE, check.names = FALSE),
            y = rewiringTable,
            by.x = joinColumn,
            by.y = joinColumn,
            all.x = TRUE,
            sort = FALSE
        )
    }

    S4Vectors::DataFrame(nodeTable, check.names = FALSE)
}

.resolveStoredOrSuppliedEdges <- function(analysisData, suppliedEdges, slotName) {
    if (!is.null(suppliedEdges)) {
        return(.flattenEdgeTables(suppliedEdges))
    }
    if (is.null(analysisData)) {
        return(NULL)
    }

    switch(
        EXPR = slotName,
        knowledgeNetworks = knowledgeNetworks(analysisData),
        correlationNetworks = correlationNetworks(analysisData),
        sparseWithinOmicCorrelations = sparseWithinOmicCorrelations(analysisData),
        sparseCrossOmicCorrelations = sparseCrossOmicCorrelations(analysisData),
        differentialCorrelationResults = differentialCorrelationResults(analysisData),
        integratedNetworks = integratedNetworks(analysisData),
        NULL
    )
}

#' Create an Integrated Multi-Omic Network
#'
#' Combine knowledge edges, dense correlation edges, sparse correlation
#' edges, and differential correlation edges into one standardized edge
#' table.
#'
#' @param analysisData Optional `MultiAssayExperiment` used to retrieve
#'   stored CorNetto results and optionally store the integrated network.
#' @param knowledgeNetwork Optional knowledge-network edge table or list.
#' @param correlationNetwork Optional dense correlation edge table or
#'   list.
#' @param sparseWithinOmicCorrelation Optional sparse within-omic edge
#'   table or list.
#' @param sparseCrossOmicCorrelation Optional sparse cross-omic edge
#'   table or list.
#' @param differentialCorrelationNetwork Optional differential
#'   correlation edge table or list.
#' @param includeReverseEdges Whether to duplicate undirected edges to
#'   support directed graph traversal.
#' @param resultName Name used when storing the integrated network.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
createIntegratedNetwork <- function(
    analysisData = NULL,
    knowledgeNetwork = NULL,
    correlationNetwork = NULL,
    sparseWithinOmicCorrelation = NULL,
    sparseCrossOmicCorrelation = NULL,
    differentialCorrelationNetwork = NULL,
    includeReverseEdges = TRUE,
    resultName = "integratedNetwork",
    storeResult = FALSE
) {
    if (!is.null(analysisData)) {
        analysisData <- validateAnalysisData(analysisData)
    }

    edgeInputs <- c(
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, knowledgeNetwork, "knowledgeNetworks")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, correlationNetwork, "correlationNetworks")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, sparseWithinOmicCorrelation, "sparseWithinOmicCorrelations")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, sparseCrossOmicCorrelation, "sparseCrossOmicCorrelations")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, differentialCorrelationNetwork, "differentialCorrelationResults"))
    )
    edgeInputs <- Filter(Negate(is.null), edgeInputs)

    if (!length(edgeInputs)) {
        stop("No edge tables were supplied for integration.", call. = FALSE)
    }

    integratedEdgeTable <- do.call(rbind, lapply(edgeInputs, .coerceStandardEdgeTable))
    integratedNodeTable <- .mergeNodeTables(edgeInputs, edgeTable = integratedEdgeTable)
    integratedEdgeTable <- .coerceStandardEdgeTable(
        integratedEdgeTable,
        nodeTable = integratedNodeTable
    )

    if (isTRUE(includeReverseEdges)) {
        integratedEdgeTable <- .duplicateUndirectedEdges(integratedEdgeTable)
    }

    if (!storeResult) {
        return(integratedEdgeTable)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "integratedNetworks",
        resultObject = integratedEdgeTable,
        resultName = resultName
    )
}

#' Filter a Network by Node Membership
#'
#' Retain only edges touching selected nodes or only edges where both
#' endpoints belong to a selected node set.
#'
#' @param networkEdgeTable A standardized edge table.
#' @param nodes Node identifiers to retain.
#' @param mode Either `"either"` or `"both"`.
#'
#' @return A standardized edge `DataFrame`.
#' @export
filterNetworkByNodes <- function(networkEdgeTable, nodes, mode = c("either", "both")) {
    mode <- match.arg(mode)
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)

    nodes <- as.character(nodes)
    inputNodeTable <- .getStoredNodeTable(networkEdgeTable, fallbackToEdges = TRUE)
    selectedNodeTable <- inputNodeTable[inputNodeTable$nodeIdentifier %in% nodes, , drop = FALSE]
    edgeData <- as.data.frame(networkEdgeTable, stringsAsFactors = FALSE, check.names = FALSE)

    keepIndex <- if (identical(mode, "either")) {
        edgeData$fromFeatureIdentifier %in% nodes | edgeData$toFeatureIdentifier %in% nodes
    } else {
        edgeData$fromFeatureIdentifier %in% nodes & edgeData$toFeatureIdentifier %in% nodes
    }

    filteredEdgeTable <- .coerceStandardEdgeTable(edgeData[keepIndex, , drop = FALSE])
    filteredNodeTable <- .mergeNodeTables(selectedNodeTable, edgeTable = filteredEdgeTable)

    .coerceStandardEdgeTable(filteredEdgeTable, nodeTable = filteredNodeTable)
}

#' Create a Focused Network Around Seed Nodes
#'
#' Build a neighborhood-induced subnetwork around a set of seed nodes.
#' Focused-network outputs retain an explicit node table in metadata so
#' isolated retained nodes remain available to downstream graph and
#' export helpers.
#'
#' @param networkEdgeTable A standardized edge table.
#' @param seedNodes Seed-node identifiers.
#' @param neighborhoodOrder Number of graph steps to expand.
#' @param mode Neighborhood mode passed to `igraph::ego()`. One of
#'   `"all"`, `"out"`, or `"in"`.
#' @param dropIsolatedNodes Whether to remove nodes retained in the
#'   focused subgraph that do not participate in any of the returned
#'   edges.
#'
#' @return A standardized edge `DataFrame`.
#' @export
createFocusedNetwork <- function(
    networkEdgeTable,
    seedNodes,
    neighborhoodOrder = 1L,
    mode = c("all", "out", "in"),
    dropIsolatedNodes = FALSE
) {
    mode <- match.arg(mode)
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)
    .assertScalarLogical(dropIsolatedNodes, "dropIsolatedNodes")

    inputNodeTable <- .getStoredNodeTable(networkEdgeTable, fallbackToEdges = TRUE)
    edgeData <- as.data.frame(networkEdgeTable, stringsAsFactors = FALSE, check.names = FALSE)
    edgeData$fromNodeKey <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    edgeData$toNodeKey <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)
    graphVertices <- as.data.frame(inputNodeTable, stringsAsFactors = FALSE, check.names = FALSE)
    graphVertices$name <- graphVertices$nodeKey

    graphObject <- igraph::graph_from_data_frame(
        d = data.frame(from = edgeData$fromNodeKey, to = edgeData$toNodeKey, stringsAsFactors = FALSE),
        directed = any(as.logical(edgeData$isDirected), na.rm = TRUE),
        vertices = graphVertices
    )

    seedNodeKeys <- unique(
        inputNodeTable$nodeKey[inputNodeTable$nodeIdentifier %in% seedNodes]
    )
    if (!length(seedNodeKeys)) {
        return(.coerceStandardEdgeTable(
            .emptyStandardEdgeTable(),
            nodeTable = .emptyStandardNodeTable()
        ))
    }

    egoNodes <- igraph::ego(
        graph = graphObject,
        order = neighborhoodOrder,
        nodes = seedNodeKeys,
        mode = mode
    )
    expandedKeys <- unique(unlist(lapply(egoNodes, igraph::as_ids), use.names = FALSE))
    keepIndex <- edgeData$fromNodeKey %in% expandedKeys & edgeData$toNodeKey %in% expandedKeys
    focusedNodeTable <- inputNodeTable[inputNodeTable$nodeKey %in% expandedKeys, , drop = FALSE]
    focusedEdgeTable <- .coerceStandardEdgeTable(
        edgeData[keepIndex, .standardEdgeColumns(), drop = FALSE]
    )

    if (isTRUE(dropIsolatedNodes)) {
        endpointNodeTable <- .getStoredNodeTable(focusedEdgeTable, fallbackToEdges = TRUE)
        focusedNodeTable <- focusedNodeTable[
            focusedNodeTable$nodeKey %in% endpointNodeTable$nodeKey,
            ,
            drop = FALSE
        ]
    }

    .coerceStandardEdgeTable(focusedEdgeTable, nodeTable = focusedNodeTable)
}
