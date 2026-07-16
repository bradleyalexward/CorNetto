.deduplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    edgeTable <- .asPlainDataFrame(edgeTable)
    keyOne <- .nodeKey(edgeTable$fromAssayName, edgeTable$fromFeatureIdentifier)
    keyTwo <- .nodeKey(edgeTable$toAssayName, edgeTable$toFeatureIdentifier)
    edgeIsDirected <- as.logical(edgeTable$isDirected)
    edgeIsDirected[is.na(edgeIsDirected)] <- FALSE

    # A->B and B->A are the same undirected edge, so order the endpoints before
    # building the key. Directed edges keep the order they were given.
    pairKey <- ifelse(
        edgeIsDirected,
        paste(keyOne, keyTwo, sep = "||"),
        paste(pmin(keyOne, keyTwo), pmax(keyOne, keyTwo), sep = "||")
    )
    pairKey <- paste(
        pairKey,
        edgeTable$edgeDirection,
        edgeTable$edgeType,
        sep = "##"
    )

    edgeTable <- edgeTable[!duplicated(pairKey), , drop = FALSE]
    .coerceStandardEdgeTable(edgeTable, nodeTable = nodeTable)
}

.selectSparseCandidateEdges <- function(
    analysisData,
    knowledgeNetwork,
    assayNames = NULL,
    correlationScope = c("all", "withinOmic", "crossOmic"),
    assayFilterMode = c("both", "either", "from", "to")
) {
    correlationScope <- match.arg(correlationScope)
    assayFilterMode <- match.arg(assayFilterMode)
    candidateEdges <- validateKnowledgeNetwork(
        knowledgeNetwork,
        storeResult = FALSE,
        quiet = TRUE
    )
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))

    candidateEdges <- candidateEdges[
        candidateEdges$fromAssayName %in% measuredAssays &
            candidateEdges$toAssayName %in% measuredAssays,
        ,
        drop = FALSE
    ]

    if (!is.null(assayNames)) {
        assayNames <- as.character(assayNames)
        keepAssay <- switch(
            EXPR = assayFilterMode,
            both = candidateEdges$fromAssayName %in% assayNames &
                candidateEdges$toAssayName %in% assayNames,
            either = candidateEdges$fromAssayName %in% assayNames |
                candidateEdges$toAssayName %in% assayNames,
            from = candidateEdges$fromAssayName %in% assayNames,
            to = candidateEdges$toAssayName %in% assayNames
        )
        candidateEdges <- candidateEdges[keepAssay, , drop = FALSE]
    }

    if (identical(correlationScope, "withinOmic")) {
        candidateEdges <- candidateEdges[
            candidateEdges$fromAssayName == candidateEdges$toAssayName,
            ,
            drop = FALSE
        ]
    } else if (identical(correlationScope, "crossOmic")) {
        candidateEdges <- candidateEdges[
            candidateEdges$fromAssayName != candidateEdges$toAssayName,
            ,
            drop = FALSE
        ]
    }

    .deduplicateUndirectedEdges(candidateEdges)
}

.runCorrelationTest <- function(x, y, correlationMethod, minimumSampleCount = 3L) {
    completeIndex <- stats::complete.cases(x, y)
    if (sum(completeIndex) < minimumSampleCount) {
        return(NULL)
    }

    x <- x[completeIndex]
    y <- y[completeIndex]
    # cor.test() has no answer for a constant vector. That is the only failure
    # we expect here, so check for it and let anything else surface.
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
        return(NULL)
    }

    testArguments <- list(x = x, y = y, method = correlationMethod)
    if (!identical(correlationMethod, "pearson")) {
        testArguments$exact <- FALSE
    }
    testResult <- suppressWarnings(do.call(stats::cor.test, testArguments))

    list(
        correlationValue = unname(testResult$estimate[[1L]]),
        pValue = unname(testResult$p.value),
        sampleCount = sum(completeIndex)
    )
}
