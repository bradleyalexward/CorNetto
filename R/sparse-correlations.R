.deduplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    edgeTable <- .asPlainDataFrame(edgeTable)
    keyOne <- paste(edgeTable$fromAssayName, edgeTable$fromFeatureIdentifier, sep = "::")
    keyTwo <- paste(edgeTable$toAssayName, edgeTable$toFeatureIdentifier, sep = "::")
    edgeIsDirected <- as.logical(edgeTable$isDirected)
    edgeIsDirected[is.na(edgeIsDirected)] <- FALSE
    pairKey <- ifelse(
        edgeIsDirected,
        paste(keyOne, keyTwo, sep = "||"),
        ifelse(keyOne <= keyTwo,
        paste(keyOne, keyTwo, sep = "||"),
        paste(keyTwo, keyOne, sep = "||")
        )
    )
    pairKey <- paste(
        pairKey,
        edgeTable$edgeDirection,
        edgeTable$edgeType,
        sep = "##"
    )
    edgeTable <- edgeTable[!duplicated(pairKey), , drop = FALSE]
    .coerceStandardEdgeTable(
        S4Vectors::DataFrame(edgeTable, check.names = FALSE),
        nodeTable = nodeTable
    )
}

.sparseCorrelationScopeForEdge <- function(fromAssayName, toAssayName) {
    if (identical(fromAssayName, toAssayName)) {
        return("sparseWithinOmic")
    }

    "sparseCrossOmic"
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

.resolveSparseFeatureName <- function(
    analysisData,
    assayName,
    featureIdentifier,
    suppliedFeatureName,
    featureNameColumn
) {
    if (!is.na(suppliedFeatureName) && nzchar(suppliedFeatureName)) {
        return(suppliedFeatureName)
    }

    featureAnnotations <- .extractFeatureAnnotations(
        analysisData = analysisData,
        assayName = assayName,
        featureNameColumn = featureNameColumn
    )
    matchedName <- featureAnnotations$featureName[
        match(featureIdentifier, featureAnnotations$featureIdentifier)
    ]
    if (!is.na(matchedName) && nzchar(matchedName)) {
        return(matchedName)
    }

    featureIdentifier
}

.runCorrelationTest <- function(x, y, correlationMethod, minimumSampleCount = 3L) {
    completeIndex <- stats::complete.cases(x, y)
    if (sum(completeIndex) < minimumSampleCount) {
        return(NULL)
    }

    testArguments <- list(
        x = x[completeIndex],
        y = y[completeIndex],
        method = correlationMethod
    )
    if (!identical(correlationMethod, "pearson")) {
        testArguments$exact <- FALSE
    }

    testResult <- tryCatch(
        suppressWarnings(do.call(stats::cor.test, testArguments)),
        error = function(...) NULL
    )
    if (is.null(testResult)) {
        return(NULL)
    }

    list(
        correlationValue = unname(testResult$estimate[[1L]]),
        pValue = unname(testResult$p.value),
        sampleCount = sum(completeIndex)
    )
}

.computeSparsePairCorrelations <- function(
    analysisData,
    candidateEdges,
    sampleIds,
    correlationMethod,
    minimumAbsoluteCorrelation,
    adjustedPValueThreshold,
    pAdjustMethod,
    featureNameColumn,
    groupName
) {
    candidateEdges <- .deduplicateUndirectedEdges(candidateEdges)
    if (!nrow(candidateEdges)) {
        return(.emptyStandardEdgeTable())
    }

    candidateEdges <- .asPlainDataFrame(candidateEdges)
    resultList <- vector("list", nrow(candidateEdges))

    for (edgeIndex in seq_len(nrow(candidateEdges))) {
        edgeRow <- candidateEdges[edgeIndex, , drop = FALSE]
        fromMatrix <- .extractAssayMatrix(analysisData, edgeRow$fromAssayName)
        toMatrix <- .extractAssayMatrix(analysisData, edgeRow$toAssayName)

        if (!edgeRow$fromFeatureIdentifier %in% rownames(fromMatrix)) {
            next
        }
        if (!edgeRow$toFeatureIdentifier %in% rownames(toMatrix)) {
            next
        }

        sharedSamples <- Reduce(
            intersect,
            list(sampleIds, colnames(fromMatrix), colnames(toMatrix))
        )
        if (length(sharedSamples) < 3L) {
            next
        }

        x <- as.numeric(fromMatrix[edgeRow$fromFeatureIdentifier, sharedSamples, drop = TRUE])
        y <- as.numeric(toMatrix[edgeRow$toFeatureIdentifier, sharedSamples, drop = TRUE])
        testResult <- .runCorrelationTest(
            x = x,
            y = y,
            correlationMethod = correlationMethod,
            minimumSampleCount = 3L
        )
        if (is.null(testResult)) {
            next
        }

        correlationValue <- testResult$correlationValue
        pValue <- testResult$pValue

        resultList[[edgeIndex]] <- data.frame(
            fromFeatureIdentifier = edgeRow$fromFeatureIdentifier,
            toFeatureIdentifier = edgeRow$toFeatureIdentifier,
            fromFeatureName = .resolveSparseFeatureName(
                analysisData = analysisData,
                assayName = edgeRow$fromAssayName,
                featureIdentifier = edgeRow$fromFeatureIdentifier,
                suppliedFeatureName = edgeRow$fromFeatureName,
                featureNameColumn = featureNameColumn
            ),
            toFeatureName = .resolveSparseFeatureName(
                analysisData = analysisData,
                assayName = edgeRow$toAssayName,
                featureIdentifier = edgeRow$toFeatureIdentifier,
                suppliedFeatureName = edgeRow$toFeatureName,
                featureNameColumn = featureNameColumn
            ),
            fromAssayName = edgeRow$fromAssayName,
            toAssayName = edgeRow$toAssayName,
            edgeType = "correlation",
            edgeDirection = ifelse(correlationValue >= 0, "positive", "negative"),
            sourceType = "sparseCorrelation",
            correlationScope = .sparseCorrelationScopeForEdge(
                edgeRow$fromAssayName,
                edgeRow$toAssayName
            ),
            correlationMethod = correlationMethod,
            knowledgeSource = edgeRow$knowledgeSource,
            groupName = groupName,
            comparisonName = NA_character_,
            correlationValue = correlationValue,
            group1CorrelationValue = NA_real_,
            group2CorrelationValue = NA_real_,
            pValue = pValue,
            adjustedPValue = NA_real_,
            group1PValue = NA_real_,
            group2PValue = NA_real_,
            group1AdjustedPValue = NA_real_,
            group2AdjustedPValue = NA_real_,
            zScoreDifference = NA_real_,
            sampleCount = testResult$sampleCount,
            group1SampleCount = NA_real_,
            group2SampleCount = NA_real_,
            edgeWeight = correlationValue,
            evidenceScore = edgeRow$evidenceScore,
            isDirected = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    resultList <- Filter(Negate(is.null), resultList)
    if (!length(resultList)) {
        return(.emptyStandardEdgeTable())
    }

    resultTable <- do.call(rbind, resultList)
    resultTable$adjustedPValue <- .adjustPValuesWithMissing(resultTable$pValue, pAdjustMethod)

    keepIndex <- !is.na(resultTable$correlationValue)
    if (!is.null(minimumAbsoluteCorrelation)) {
        keepIndex <- keepIndex & abs(resultTable$correlationValue) >= minimumAbsoluteCorrelation
    }
    if (!is.null(adjustedPValueThreshold)) {
        keepIndex <- keepIndex & resultTable$adjustedPValue <= adjustedPValueThreshold
    }

    .coerceStandardEdgeTable(resultTable[keepIndex, , drop = FALSE])
}

.resolveSparseGroupName <- function(groupColumn, groupLevel, sampleIds) {
    if (!is.null(groupLevel)) {
        return(groupLevel)
    }
    if (!is.null(groupColumn)) {
        return(groupColumn)
    }
    .makeResultName("customSamples", length(sampleIds))
}
