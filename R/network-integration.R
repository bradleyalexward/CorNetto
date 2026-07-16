.nodeKey <- function(assayName, featureIdentifier) {
    assayName <- as.character(assayName)
    featureIdentifier <- as.character(featureIdentifier)
    hasSeparator <- (!is.na(assayName) & grepl("::", assayName, fixed = TRUE)) |
        (!is.na(featureIdentifier) &
            grepl("::", featureIdentifier, fixed = TRUE))
    if (any(hasSeparator)) {
        stop(
            "Assay names and feature identifiers cannot contain `::`.",
            call. = FALSE
        )
    }
    paste(assayName, featureIdentifier, sep = "::")
}

## igraph cannot hold directed and undirected edges in one object: if any edge
## is directed the whole graph is. Mirroring the undirected edges first is what
## keeps them symmetric under mode = "in"/"out" traversal.
.graphIsDirected <- function(edgeData) {
    if (!nrow(edgeData)) {
        return(FALSE)
    }

    any(as.logical(edgeData$isDirected), na.rm = TRUE)
}

.duplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    edgeData <- .asPlainDataFrame(edgeTable)
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

    # Mirroring an already-mirrored table must be a no-op, otherwise composing
    # combineNetworks() with createNetworkGraph() gives four copies of every
    # undirected edge. The mirror of a mirror is the original row, so dropping
    # reverses that are already present makes this idempotent.
    if (nrow(duplicateTable)) {
        alreadyPresent <- .edgeRowKey(duplicateTable) %in% .edgeRowKey(edgeData)
        duplicateTable <- duplicateTable[!alreadyPresent, , drop = FALSE]
    }
    if (!nrow(duplicateTable)) {
        return(edgeTable)
    }

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
        rewiringTable <- .asPlainDataFrame(rewiringTable)
        if (!any(c("nodeKey", "nodeIdentifier") %in% names(rewiringTable))) {
            stop("`rewiringTable` must contain `nodeKey` or `nodeIdentifier`.", call. = FALSE)
        }
        joinColumn <- if ("nodeKey" %in% names(rewiringTable)) {
            "nodeKey"
        } else {
            if (anyDuplicated(nodeTable$nodeIdentifier)) {
                stop(
                    "`rewiringTable` must contain `nodeKey` when feature identifiers occur in multiple assays.",
                    call. = FALSE
                )
            }
            "nodeIdentifier"
        }
        if (anyDuplicated(rewiringTable[[joinColumn]])) {
            stop(
                "`rewiringTable` must contain at most one row per `",
                joinColumn,
                "`.",
                call. = FALSE
            )
        }
        rewiringTable <- rewiringTable[
            ,
            setdiff(names(rewiringTable), setdiff(.standardNodeColumns(), joinColumn)),
            drop = FALSE
        ]
        nodeTable <- merge(
            x = .asPlainDataFrame(nodeTable),
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
        correlationResults = correlationResults(analysisData),
        differentialCorrelationResults = differentialCorrelationResults(analysisData),
        differentialCorrelationNetworks = differentialCorrelationNetworks(analysisData),
        integratedNetworks = integratedNetworks(analysisData),
        NULL
    )
}

#' Combine Network Edge Tables
#'
#' Combine knowledge edges, correlation edges, differential-correlation
#' edges, or other standardized CorNetto edge tables into one network.
#'
#' @param analysisData Optional `MultiAssayExperiment` used to retrieve
#'   stored CorNetto results and optionally store the combined network.
#' @param knowledgeNetwork Optional knowledge-network edge table or list
#'   of edge tables.
#' @param correlationNetwork Optional correlation edge table or list of
#'   edge tables.
#' @param differentialCorrelationNetwork Optional differential
#'   correlation edge table or list of edge tables.
#' @param networkList Optional additional edge table or list of edge
#'   tables.
#' @param includeReverseEdges Whether to duplicate undirected edges to
#'   support directed graph traversal. See Details before enabling it for a
#'   table you intend to score.
#' @param removeDuplicates Whether to remove exactly duplicated edge
#'   rows after standardization.
#' @param resultName Name used when storing the combined network.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @details `includeReverseEdges = TRUE` appends a reversed copy of every
#'   undirected edge. That is what you want for traversal and for export to
#'   tools that expect an explicit adjacency list, and it is the default for
#'   that reason.
#'
#'   It is not what you want for anything that treats the table as a set of
#'   edges. [calculateRewiringScores()] rejects duplicate or mirrored edges
#'   because they would double every node's degree and multiply
#'   `rawRewiringScore` by `sqrt(2)`. An undirected [createNetworkGraph()] built
#'   from a mirrored table gains parallel edges.
#'   Score first, then combine with reverse edges for export, or pass
#'   `includeReverseEdges = FALSE` and let [createNetworkGraph()] handle
#'   direction.
#' @examples
#' analysisData <- exampleAnalysisData()
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' differentialResults <- testDifferentialCorrelation(
#'     analysisData = analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
#'     storeResult = FALSE
#' )
#' differentialNetwork <- createDifferentialCorrelationNetwork(
#'     differentialResults,
#'     minimumAbsoluteCorrelation = 0
#' )
#' combinedNetwork <- combineNetworks(
#'     knowledgeNetwork = knowledgeNetwork,
#'     differentialCorrelationNetwork = differentialNetwork
#' )
#' combinedNetwork
#' @export
combineNetworks <- function(
    analysisData = NULL,
    knowledgeNetwork = NULL,
    correlationNetwork = NULL,
    differentialCorrelationNetwork = NULL,
    networkList = NULL,
    includeReverseEdges = TRUE,
    removeDuplicates = TRUE,
    resultName = "combinedNetwork",
    storeResult = FALSE
) {
    if (!is.null(analysisData)) {
        analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    }
    .assertScalarLogical(includeReverseEdges, "includeReverseEdges")
    .assertScalarLogical(removeDuplicates, "removeDuplicates")

    edgeInputs <- c(
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, knowledgeNetwork, "knowledgeNetworks")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, correlationNetwork, "correlationResults")),
        .flattenEdgeTables(.resolveStoredOrSuppliedEdges(analysisData, differentialCorrelationNetwork, "differentialCorrelationNetworks")),
        .flattenEdgeTables(networkList)
    )
    edgeInputs <- Filter(Negate(is.null), edgeInputs)

    if (!length(edgeInputs)) {
        stop("No edge tables were supplied to combine.", call. = FALSE)
    }

    integratedEdgeTable <- do.call(rbind, lapply(edgeInputs, .coerceStandardEdgeTable))
    integratedNodeTable <- .mergeNodeTables(edgeInputs, edgeTable = integratedEdgeTable)
    if (isTRUE(removeDuplicates) && nrow(integratedEdgeTable)) {
        integratedEdgeTable <- unique(.asPlainDataFrame(integratedEdgeTable))
    }
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
#' @param nodes Node identifiers or node keys to retain.
#' @param mode Either `"either"` or `"both"`.
#' @param nodeIdType Whether `nodes` are assay-local feature
#'   identifiers (`"identifier"`) or assay-qualified node keys
#'   (`"key"`).
#'
#' @return A standardized edge `DataFrame`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' filterNetworkByNodes(knowledgeNetwork, nodes = c("P1", "P2"))
#' @export
filterNetworkByNodes <- function(
    networkEdgeTable,
    nodes,
    mode = c("either", "both"),
    nodeIdType = c("identifier", "key")
) {
    mode <- match.arg(mode)
    nodeIdType <- match.arg(nodeIdType)
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)

    nodes <- as.character(nodes)
    inputNodeTable <- .getStoredNodeTable(networkEdgeTable, fallbackToEdges = TRUE)
    selectedColumn <- if (identical(nodeIdType, "key")) "nodeKey" else "nodeIdentifier"
    selectedNodeTable <- inputNodeTable[inputNodeTable[[selectedColumn]] %in% nodes, , drop = FALSE]
    edgeData <- .asPlainDataFrame(networkEdgeTable)
    fromNodeKey <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    toNodeKey <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)

    keepIndex <- if (identical(mode, "either")) {
        if (identical(nodeIdType, "key")) {
            fromNodeKey %in% nodes | toNodeKey %in% nodes
        } else {
            edgeData$fromFeatureIdentifier %in% nodes | edgeData$toFeatureIdentifier %in% nodes
        }
    } else {
        if (identical(nodeIdType, "key")) {
            fromNodeKey %in% nodes & toNodeKey %in% nodes
        } else {
            edgeData$fromFeatureIdentifier %in% nodes & edgeData$toFeatureIdentifier %in% nodes
        }
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
#' @param seedNodes Seed-node identifiers or node keys.
#' @param neighborhoodOrder Number of graph steps to expand.
#' @param mode Neighborhood mode passed to `igraph::ego()`. One of
#'   `"all"`, `"out"`, or `"in"`.
#' @param nodeIdType Whether `seedNodes` are assay-local feature
#'   identifiers (`"identifier"`) or assay-qualified node keys
#'   (`"key"`).
#' @param dropIsolatedNodes Whether to remove nodes retained in the
#'   focused subgraph that do not participate in any of the returned
#'   edges.
#'
#' @return A standardized edge `DataFrame`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' focusedNetwork <- createFocusedNetwork(
#'     knowledgeNetwork,
#'     seedNodes = exampleSeedNodes()
#' )
#' focusedNetwork
#' @export
createFocusedNetwork <- function(
    networkEdgeTable,
    seedNodes,
    neighborhoodOrder = 1L,
    mode = c("all", "out", "in"),
    nodeIdType = c("identifier", "key"),
    dropIsolatedNodes = FALSE
) {
    mode <- match.arg(mode)
    nodeIdType <- match.arg(nodeIdType)
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)
    .assertScalarLogical(dropIsolatedNodes, "dropIsolatedNodes")
    neighborhoodOrder <- .assertWholeNumber(
        neighborhoodOrder,
        "neighborhoodOrder",
        minimum = 0L
    )

    inputNodeTable <- .getStoredNodeTable(networkEdgeTable, fallbackToEdges = TRUE)
    edgeData <- .asPlainDataFrame(networkEdgeTable)
    edgeData$fromNodeKey <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    edgeData$toNodeKey <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)
    graphVertices <- .asPlainDataFrame(inputNodeTable)
    graphVertices$name <- graphVertices$nodeKey

    # Traverse on a mirrored copy when the graph is directed, so an undirected
    # edge is reachable from both ends. The returned edges come from edgeData,
    # so this never leaks duplicates into the result.
    graphDirected <- .graphIsDirected(edgeData)
    traversalEdges <- if (graphDirected) {
        .asPlainDataFrame(.duplicateUndirectedEdges(networkEdgeTable))
    } else {
        edgeData
    }

    graphObject <- igraph::graph_from_data_frame(
        d = data.frame(
            from = .nodeKey(traversalEdges$fromAssayName, traversalEdges$fromFeatureIdentifier),
            to = .nodeKey(traversalEdges$toAssayName, traversalEdges$toFeatureIdentifier),
            stringsAsFactors = FALSE
        ),
        directed = graphDirected,
        vertices = graphVertices
    )

    seedColumn <- if (identical(nodeIdType, "key")) "nodeKey" else "nodeIdentifier"
    seedNodeKeys <- unique(inputNodeTable$nodeKey[inputNodeTable[[seedColumn]] %in% seedNodes])
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
