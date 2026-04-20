.deduplicateUndirectedEdges <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    nodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    edgeTable <- as.data.frame(edgeTable, stringsAsFactors = FALSE, check.names = FALSE)
    keyOne <- paste(edgeTable$fromAssayName, edgeTable$fromFeatureIdentifier, sep = "::")
    keyTwo <- paste(edgeTable$toAssayName, edgeTable$toFeatureIdentifier, sep = "::")
    pairKey <- ifelse(keyOne <= keyTwo, paste(keyOne, keyTwo, sep = "||"), paste(keyTwo, keyOne, sep = "||"))
    edgeTable <- edgeTable[!duplicated(pairKey), , drop = FALSE]
    .coerceStandardEdgeTable(
        S4Vectors::DataFrame(edgeTable, check.names = FALSE),
        nodeTable = nodeTable
    )
}

.computeSparsePairCorrelations <- function(
    analysisData,
    candidateEdges,
    sampleIds,
    correlationMethod,
    minimumAbsoluteCorrelation,
    adjustedPValueThreshold,
    featureNameColumn,
    correlationScope,
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
        completeIndex <- stats::complete.cases(x, y)
        if (sum(completeIndex) < 3L) {
            next
        }

        testResult <- tryCatch(
            suppressWarnings(
                stats::cor.test(
                    x = x[completeIndex],
                    y = y[completeIndex],
                    method = correlationMethod,
                    exact = FALSE
                )
            ),
            error = function(...) NULL
        )
        if (is.null(testResult)) {
            next
        }

        correlationValue <- unname(testResult$estimate[[1L]])
        pValue <- unname(testResult$p.value)

        resultList[[edgeIndex]] <- data.frame(
            fromFeatureIdentifier = edgeRow$fromFeatureIdentifier,
            toFeatureIdentifier = edgeRow$toFeatureIdentifier,
            fromFeatureName = if (!is.na(edgeRow$fromFeatureName)) edgeRow$fromFeatureName else edgeRow$fromFeatureIdentifier,
            toFeatureName = if (!is.na(edgeRow$toFeatureName)) edgeRow$toFeatureName else edgeRow$toFeatureIdentifier,
            fromAssayName = edgeRow$fromAssayName,
            toAssayName = edgeRow$toAssayName,
            edgeType = "correlation",
            edgeDirection = ifelse(correlationValue >= 0, "positive", "negative"),
            sourceType = "sparseCorrelation",
            correlationScope = correlationScope,
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
    resultTable$adjustedPValue <- stats::p.adjust(resultTable$pValue, method = "fdr")

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

#' Create Sparse Within-Omic Correlations from Prior Edges
#'
#' Compute within-omic correlations only for feature pairs supported by a
#' prior-knowledge network.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param knowledgeNetwork A standardized knowledge network or any
#'   table accepted by `validateKnowledgeNetwork()`.
#' @param assayNames Optional assay names to retain. When `NULL`, all
#'   within-omic assay pairs in the prior network are considered.
#' @param groupColumn Sample metadata column used to define the group.
#' @param groupLevel Single group label to analyse.
#' @param sampleIds Optional explicit sample identifiers.
#' @param correlationMethod Correlation method.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation
#'   retained.
#' @param adjustedPValueThreshold Maximum adjusted p-value retained.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional storage name.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
createSparseWithinOmicCorrelations <- function(
    analysisData,
    knowledgeNetwork,
    assayNames = NULL,
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    featureNameColumn = "featureName",
    resultName = NULL,
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationMethod <- match.arg(correlationMethod)
    candidateEdges <- validateKnowledgeNetwork(knowledgeNetwork, storeResult = FALSE)
    candidateEdges <- candidateEdges[
        candidateEdges$fromAssayName == candidateEdges$toAssayName,
        ,
        drop = FALSE
    ]

    if (!is.null(assayNames)) {
        candidateEdges <- candidateEdges[
            candidateEdges$fromAssayName %in% assayNames,
            ,
            drop = FALSE
        ]
    }

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
        featureNameColumn = featureNameColumn,
        correlationScope = "sparseWithinOmic",
        groupName = groupName
    )

    if (!storeResult) {
        return(edgeTable)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName("sparseWithinOmic", groupName, correlationMethod)
    }

    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "sparseWithinOmicCorrelations",
        resultObject = edgeTable,
        resultName = resultName
    )
}

#' Create Sparse Cross-Omic Correlations from Prior Edges
#'
#' Compute cross-omic correlations only for feature pairs supported by a
#' prior-knowledge network.
#'
#' @inheritParams createSparseWithinOmicCorrelations
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
createSparseCrossOmicCorrelations <- function(
    analysisData,
    knowledgeNetwork,
    assayNames = NULL,
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    featureNameColumn = "featureName",
    resultName = NULL,
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationMethod <- match.arg(correlationMethod)
    candidateEdges <- validateKnowledgeNetwork(knowledgeNetwork, storeResult = FALSE)
    candidateEdges <- candidateEdges[
        candidateEdges$fromAssayName != candidateEdges$toAssayName,
        ,
        drop = FALSE
    ]

    if (!is.null(assayNames)) {
        candidateEdges <- candidateEdges[
            candidateEdges$fromAssayName %in% assayNames |
                candidateEdges$toAssayName %in% assayNames,
            ,
            drop = FALSE
        ]
    }

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
        featureNameColumn = featureNameColumn,
        correlationScope = "sparseCrossOmic",
        groupName = groupName
    )

    if (!storeResult) {
        return(edgeTable)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName("sparseCrossOmic", groupName, correlationMethod)
    }

    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "sparseCrossOmicCorrelations",
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
        undirectedKey <- ifelse(
            combinedData$isDirected,
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
