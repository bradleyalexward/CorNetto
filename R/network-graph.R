#' Create an igraph Object from a CorNetto Edge Table
#'
#' Stored node metadata is used when present so isolated nodes can be
#' preserved across graph creation and export.
#'
#' @param networkEdgeTable A standardized edge table.
#' @param nodeTable Optional node table. When `NULL`, the stored node
#'   table is used when present and otherwise derived from the edges.
#' @param directed Optional logical value. When `NULL`, the graph is
#'   directed when any edge is directed.
#' @param mirrorUndirectedEdges Whether to add a reverse copy of each
#'   undirected edge when the graph is directed. `igraph` cannot mix
#'   directed and undirected edges in one object, so without this an
#'   undirected edge in a directed graph is only traversable one way.
#'
#' @return An `igraph` object.
#' @details `igraph` graphs are either directed or undirected as a whole. When
#'   a table mixes both, the graph is built directed and each undirected edge
#'   is added in both directions so that `mode = "in"` and `mode = "out"`
#'   traversal treat it symmetrically. Set `mirrorUndirectedEdges = FALSE` to
#'   keep a one-to-one correspondence between rows and graph edges, at the cost
#'   of asymmetric neighbourhoods.
#'
#'   All standardized edge columns become edge attributes. Nonstandard columns
#'   in a supplied or stored node table are retained as vertex attributes.
#'
#'   Note that an edge table already carrying reverse copies, as produced by
#'   `combineNetworks(includeReverseEdges = TRUE)`, will yield parallel edges
#'   in an undirected graph. Prefer `includeReverseEdges = FALSE` and let this
#'   function handle direction.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' graph <- createNetworkGraph(knowledgeNetwork)
#' graph
#' @export
createNetworkGraph <- function(
    networkEdgeTable,
    nodeTable = NULL,
    directed = NULL,
    mirrorUndirectedEdges = TRUE
) {
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)
    .assertScalarLogical(mirrorUndirectedEdges, "mirrorUndirectedEdges")
    if (!is.null(directed)) {
        .assertScalarLogical(directed, "directed")
    }
    if (is.null(nodeTable)) {
        nodeTable <- .getStoredNodeTable(
            networkEdgeTable,
            fallbackToEdges = TRUE,
            preserveExtraColumns = TRUE
        )
    }
    nodeTable <- .coerceStandardNodeTable(
        nodeTable,
        preserveExtraColumns = TRUE
    )
    nodeTable <- .asPlainDataFrame(nodeTable)
    nodeTable$name <- nodeTable$nodeKey

    graphDirected <- directed
    if (is.null(graphDirected)) {
        graphDirected <- .graphIsDirected(.asPlainDataFrame(networkEdgeTable))
    }

    if (isTRUE(graphDirected) && isTRUE(mirrorUndirectedEdges)) {
        networkEdgeTable <- .duplicateUndirectedEdges(networkEdgeTable)
    }

    edgeData <- .asPlainDataFrame(networkEdgeTable)
    edgeData$fromNodeKey <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    edgeData$toNodeKey <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)

    if (!nrow(edgeData)) {
        if (!nrow(nodeTable)) {
            stop("`networkEdgeTable` does not contain any edges or stored nodes.", call. = FALSE)
        }

        return(
            igraph::graph_from_data_frame(
                d = data.frame(from = character(), to = character(), stringsAsFactors = FALSE),
                directed = graphDirected,
                vertices = nodeTable
            )
        )
    }

    graphEdges <- data.frame(
        from = edgeData$fromNodeKey,
        to = edgeData$toNodeKey,
        edgeData[, .standardEdgeColumns(), drop = FALSE],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    igraph::graph_from_data_frame(
        d = graphEdges,
        directed = graphDirected,
        vertices = nodeTable
    )
}

.resolveRewiringPlotLabels <- function(rewiringTable, label = NULL) {
    if (is.null(label)) {
        if ("nodeKey" %in% names(rewiringTable)) {
            return(as.character(rewiringTable$nodeKey))
        }
        return(as.character(rewiringTable$nodeIdentifier))
    }

    label <- match.arg(label, c("identifier", "name", "key"))
    labelColumn <- switch(
        EXPR = label,
        identifier = "nodeIdentifier",
        name = "nodeName",
        key = "nodeKey"
    )
    if (!labelColumn %in% names(rewiringTable)) {
        stop("Column `", labelColumn, "` was not found in `rewiringTable`.", call. = FALSE)
    }

    labels <- as.character(rewiringTable[[labelColumn]])
    if (identical(label, "name") && "nodeIdentifier" %in% names(rewiringTable)) {
        missingLabel <- is.na(labels) | !nzchar(labels)
        labels[missingLabel] <- as.character(rewiringTable$nodeIdentifier[missingLabel])
    }

    labels
}

#' Plot Rewiring Scores
#'
#' @param rewiringTable Output of `calculateRewiringScores()`.
#' @param scoreColumn Score column to plot. See
#'   [calculateRewiringScores()] for what the three rewiring scores mean.
#'   The tail-probability columns from [permuteRewiringScores()] are also
#'   supported and are ordered from smallest to largest.
#' @param topN Number of top nodes to display.
#' @param label Optional label source. When `NULL`, `nodeKey` is used
#'   when present and `nodeIdentifier` otherwise. Use `"identifier"`,
#'   `"name"`, or `"key"` to select a specific label source.
#' @param main Plot title.
#'
#' @return Invisibly returns the plotted subset.
#' @details The default is `rootMeanSquareRewiringScore` because it is
#'   comparable across nodes of different degree and is defined for every
#'   node. `degreeMatchedZScore` is only a within-bin ranking aid and is `NA`
#'   for sparsely populated degree bins, so it is a poor default for small
#'   networks. Rewiring scores are ordered decreasingly. Permutation tail
#'   probabilities are ordered increasingly because smaller values are more
#'   extreme.
#' @examples
#' rewiringTable <- S4Vectors::DataFrame(
#'     nodeKey = c("protein::P1", "protein::P2"),
#'     nodeIdentifier = c("P1", "P2"),
#'     nodeName = c("CFH", "C3"),
#'     rootMeanSquareRewiringScore = c(2, 1),
#'     check.names = FALSE
#' )
#' plotRewiringScores(rewiringTable, topN = 2)
#' @export
plotRewiringScores <- function(
    rewiringTable,
    scoreColumn = "rootMeanSquareRewiringScore",
    topN = 20L,
    label = NULL,
    main = "Top rewired nodes"
) {
    rewiringTable <- .asPlainDataFrame(rewiringTable)
    .assertScalarCharacter(scoreColumn, "scoreColumn")
    .assertScalarCharacter(main, "main")
    if (!scoreColumn %in% names(rewiringTable)) {
        stop("`scoreColumn` was not found in `rewiringTable`.", call. = FALSE)
    }
    topN <- .assertWholeNumber(topN, "topN", minimum = 1L)
    if (!any(!is.na(rewiringTable[[scoreColumn]]))) {
        stop("`scoreColumn` contains only missing values.", call. = FALSE)
    }

    # order() puts NAs last, so without dropping them first a table with fewer
    # than topN scored nodes plots blank bars.
    rewiringTable <- rewiringTable[!is.na(rewiringTable[[scoreColumn]]), , drop = FALSE]
    tailProbabilityColumn <- scoreColumn %in% c(
        "permutationTailProbability",
        "adjustedPermutationTailProbability"
    )
    rewiringTable <- rewiringTable[
        order(
            rewiringTable[[scoreColumn]],
            decreasing = !tailProbabilityColumn
        ),
        ,
        drop = FALSE
    ]
    rewiringTable <- utils::head(rewiringTable, topN)
    labels <- .resolveRewiringPlotLabels(rewiringTable, label = label)

    graphics::barplot(
        height = rewiringTable[[scoreColumn]],
        names.arg = labels,
        las = 2,
        col = "#2B8CBE",
        border = NA,
        main = main,
        ylab = scoreColumn
    )

    invisible(rewiringTable)
}

#' Prepare Cytoscape Node and Edge Tables
#'
#' @param networkEdgeTable A standardized edge table.
#' @param rewiringTable Optional rewiring table to merge into the node
#'   table. Stored node metadata is retained when present so isolated
#'   nodes can be exported.
#'
#' @return A named list with `nodes` and `edges`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' networkTables <- prepareCytoscapeTables(knowledgeNetwork)
#' names(networkTables)
#' @export
prepareCytoscapeTables <- function(networkEdgeTable, rewiringTable = NULL) {
    networkEdgeTable <- .coerceStandardEdgeTable(networkEdgeTable)
    nodeTable <- .createNodeTable(
        edgeTable = networkEdgeTable,
        rewiringTable = rewiringTable,
        nodeTable = .getStoredNodeTable(networkEdgeTable, fallbackToEdges = TRUE)
    )
    edgeData <- .asPlainDataFrame(networkEdgeTable)
    edgeData$source <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    edgeData$target <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)

    list(
        nodes = nodeTable,
        edges = S4Vectors::DataFrame(edgeData, check.names = FALSE)
    )
}

#' Write Cytoscape-Ready Network Tables
#'
#' @param networkTables Output of `prepareCytoscapeTables()`, or a named
#'   list containing `nodes` and `edges`.
#' @param directoryPath Output directory. Required, and created if it does
#'   not exist.
#' @param prefix File-name prefix.
#' @param fileFormat Either `"tsv"` or `"csv"`.
#'
#' @return Invisibly returns the written file paths.
#' @details `directoryPath` has no default so that the destination is always a
#'   deliberate choice rather than the current working directory.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' networkTables <- prepareCytoscapeTables(knowledgeNetwork)
#' outputDir <- file.path(tempdir(), "cornetto-example")
#' writtenFiles <- writeNetworkTables(networkTables, outputDir)
#' file.exists(writtenFiles)
#' @export
writeNetworkTables <- function(
    networkTables,
    directoryPath,
    prefix = "cornettoNetwork",
    fileFormat = c("tsv", "csv")
) {
    fileFormat <- match.arg(fileFormat)
    .assertScalarCharacter(directoryPath, "directoryPath")
    .assertScalarCharacter(prefix, "prefix")
    if (!all(c("nodes", "edges") %in% names(networkTables))) {
        stop("`networkTables` must contain `nodes` and `edges`.", call. = FALSE)
    }

    if (!dir.exists(directoryPath)) {
        dir.create(directoryPath, recursive = TRUE)
    }

    extension <- fileFormat
    nodePath <- file.path(directoryPath, paste0(prefix, "_nodes.", extension))
    edgePath <- file.path(directoryPath, paste0(prefix, "_edges.", extension))

    writer <- if (identical(fileFormat, "csv")) readr::write_csv else readr::write_tsv
    writer(.asPlainDataFrame(networkTables$nodes), nodePath)
    writer(.asPlainDataFrame(networkTables$edges), edgePath)

    invisible(c(nodes = nodePath, edges = edgePath))
}
