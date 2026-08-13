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
