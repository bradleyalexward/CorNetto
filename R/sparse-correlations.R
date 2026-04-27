.deduplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    edgeTable <- as.data.frame(edgeTable, stringsAsFactors = FALSE, check.names = FALSE)
    keyOne <- paste(edgeTable$fromAssayName, edgeTable$fromFeatureIdentifier, sep = "::")
    keyTwo <- paste(edgeTable$toAssayName, edgeTable$toFeatureIdentifier, sep = "::")
    pairKey <- ifelse(
        keyOne <= keyTwo,
        paste(keyOne, keyTwo, sep = "||"),
        paste(keyTwo, keyOne, sep = "||")
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
    correlationScope = c("all", "withinOmic", "crossOmic")
) {
    correlationScope <- match.arg(correlationScope)
    candidateEdges <- validateKnowledgeNetwork(knowledgeNetwork, storeResult = FALSE)
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))

    candidateEdges <- candidateEdges[
        candidateEdges$fromAssayName %in% measuredAssays &
            candidateEdges$toAssayName %in% measuredAssays,
        ,
        drop = FALSE
    ]

    if (!is.null(assayNames)) {
        assayNames <- as.character(assayNames)
        candidateEdges <- candidateEdges[
            candidateEdges$fromAssayName %in% assayNames |
                candidateEdges$toAssayName %in% assayNames,
            ,
            drop = FALSE
        ]
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

    candidateEdges <- as.data.frame(candidateEdges, stringsAsFactors = FALSE, check.names = FALSE)
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

#' Create Sparse Multi-Omic Correlations from Prior Edges
#'
#' Compute prior-guided correlations for within-omic, cross-omic, or all
#' measured feature pairs represented in a knowledge network.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param knowledgeNetwork A standardized knowledge network or any
#'   table accepted by `validateKnowledgeNetwork()`.
#' @param assayNames Optional assay names to retain. When `NULL`, all
#'   measured assays represented in the prior network are considered.
#' @param correlationScope Which prior edges to test: `"all"`,
#'   `"withinOmic"`, or `"crossOmic"`.
#' @param groupColumn Sample metadata column used to define the group.
#' @param groupLevel Single group label to analyse.
#' @param sampleIds Optional explicit sample identifiers.
#' @param correlationMethod Correlation method.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation
#'   retained.
#' @param adjustedPValueThreshold Maximum adjusted p-value retained.
#' @param pAdjustMethod Multiple-testing correction method. Use
#'   `"qvalue"` or any method accepted by `stats::p.adjust()`.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional storage name.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
createSparseMultiOmicCorrelations <- function(
    analysisData,
    knowledgeNetwork,
    assayNames = NULL,
    correlationScope = c("all", "withinOmic", "crossOmic"),
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    pAdjustMethod = "fdr",
    featureNameColumn = "featureName",
    resultName = NULL,
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationScope <- match.arg(correlationScope)
    correlationMethod <- match.arg(correlationMethod)
    candidateEdges <- .selectSparseCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = knowledgeNetwork,
        assayNames = assayNames,
        correlationScope = correlationScope
    )

    sampleIds <- .resolveGroupSamples(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevel = groupLevel,
        sampleIds = sampleIds
    )
    groupName <- .resolveSparseGroupName(groupColumn, groupLevel, sampleIds)
    edgeTable <- .computeSparsePairCorrelations(
        analysisData = analysisData,
        candidateEdges = candidateEdges,
        sampleIds = sampleIds,
        correlationMethod = correlationMethod,
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        adjustedPValueThreshold = adjustedPValueThreshold,
        pAdjustMethod = pAdjustMethod,
        featureNameColumn = featureNameColumn,
        groupName = groupName
    )

    if (!storeResult) {
        return(edgeTable)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName("sparseMultiOmic", correlationScope, groupName, correlationMethod)
    }

    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "sparseMultiOmicCorrelations",
        resultObject = edgeTable,
        resultName = resultName
    )
}

#' Combine Correlation Edge Tables
#'
#' Combine dense and sparse correlation edge tables into one harmonized
#' dynamic network layer.
#'
#' @param ... Edge tables or lists of edge tables.
#' @param removeDuplicates Whether to remove duplicated interactions.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   combined edge table.
#' @param resultName Name used when storing the combined result.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
combineCorrelationNetworks <- function(
    ...,
    removeDuplicates = TRUE,
    analysisData = NULL,
    resultName = "combinedCorrelationNetwork",
    storeResult = FALSE
) {
    edgeTables <- .flattenEdgeTables(...)
    if (!length(edgeTables)) {
        stop("At least one correlation edge table must be supplied.", call. = FALSE)
    }

    edgeTables <- lapply(edgeTables, .coerceStandardEdgeTable)
    combinedTable <- do.call(rbind, edgeTables)
    combinedNodeTable <- .mergeNodeTables(edgeTables, edgeTable = combinedTable)

    if (removeDuplicates && nrow(combinedTable)) {
        combinedData <- as.data.frame(combinedTable, stringsAsFactors = FALSE, check.names = FALSE)
        keyOne <- paste(combinedData$fromAssayName, combinedData$fromFeatureIdentifier, sep = "::")
        keyTwo <- paste(combinedData$toAssayName, combinedData$toFeatureIdentifier, sep = "::")
        edgeIsDirected <- as.logical(combinedData$isDirected)
        edgeIsDirected[is.na(edgeIsDirected)] <- FALSE
        undirectedKey <- ifelse(
            edgeIsDirected,
            paste(keyOne, keyTwo, sep = "||"),
            ifelse(keyOne <= keyTwo, paste(keyOne, keyTwo, sep = "||"), paste(keyTwo, keyOne, sep = "||"))
        )
        fullKey <- paste(
            undirectedKey,
            combinedData$edgeType,
            combinedData$sourceType,
            combinedData$correlationScope,
            combinedData$correlationMethod,
            combinedData$groupName,
            combinedData$comparisonName,
            sep = "##"
        )
        combinedData <- combinedData[!duplicated(fullKey), , drop = FALSE]
        combinedTable <- .coerceStandardEdgeTable(combinedData, nodeTable = combinedNodeTable)
    } else {
        combinedTable <- .coerceStandardEdgeTable(combinedTable, nodeTable = combinedNodeTable)
    }

    if (!storeResult) {
        return(combinedTable)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "correlationNetworks",
        resultObject = combinedTable,
        resultName = resultName
    )
}
