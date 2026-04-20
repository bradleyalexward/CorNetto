.computePearsonCorrelations <- function(assayMatrix) {
    if (nrow(assayMatrix) < 2L) {
        return(list(
            correlationMatrix = matrix(numeric(), nrow = 0L, ncol = 0L),
            pValueMatrix = matrix(numeric(), nrow = 0L, ncol = 0L)
        ))
    }

    correlationMatrix <- stats::cor(
        x = t(assayMatrix),
        use = "pairwise.complete.obs",
        method = "pearson"
    )
    sampleCountMatrix <- tcrossprod(!is.na(assayMatrix))
    denominator <- pmax(1e-12, 1 - (correlationMatrix ^ 2))
    statisticMatrix <- correlationMatrix * sqrt(pmax(sampleCountMatrix - 2, 0) / denominator)
    pValueMatrix <- 2 * stats::pt(
        q = abs(statisticMatrix),
        df = pmax(sampleCountMatrix - 2, 1),
        lower.tail = FALSE
    )

    pValueMatrix[sampleCountMatrix < 3] <- NA_real_
    diag(pValueMatrix) <- 0
    diag(correlationMatrix) <- 1

    list(
        correlationMatrix = correlationMatrix,
        pValueMatrix = pValueMatrix
    )
}

.computeCorrelationByPairwiseTests <- function(assayMatrix, correlationMethod) {
    featureCount <- nrow(assayMatrix)
    featureNames <- rownames(assayMatrix)
    correlationMatrix <- matrix(
        NA_real_,
        nrow = featureCount,
        ncol = featureCount,
        dimnames = list(featureNames, featureNames)
    )
    pValueMatrix <- matrix(
        NA_real_,
        nrow = featureCount,
        ncol = featureCount,
        dimnames = list(featureNames, featureNames)
    )

    diag(correlationMatrix) <- 1
    diag(pValueMatrix) <- 0

    for (rowIndex in seq_len(featureCount - 1L)) {
        for (columnIndex in seq.int(rowIndex + 1L, featureCount)) {
            x <- as.numeric(assayMatrix[rowIndex, ])
            y <- as.numeric(assayMatrix[columnIndex, ])
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

            estimate <- unname(testResult$estimate[[1L]])
            pValue <- unname(testResult$p.value)
            correlationMatrix[rowIndex, columnIndex] <- estimate
            correlationMatrix[columnIndex, rowIndex] <- estimate
            pValueMatrix[rowIndex, columnIndex] <- pValue
            pValueMatrix[columnIndex, rowIndex] <- pValue
        }
    }

    list(
        correlationMatrix = correlationMatrix,
        pValueMatrix = pValueMatrix
    )
}

.computeSpearmanCorrelations <- function(assayMatrix) {
    .computeCorrelationByPairwiseTests(assayMatrix = assayMatrix, correlationMethod = "spearman")
}

.computeKendallCorrelations <- function(assayMatrix) {
    .computeCorrelationByPairwiseTests(assayMatrix = assayMatrix, correlationMethod = "kendall")
}

.dispatchCorrelationMethod <- function(assayMatrix, correlationMethod) {
    correlationMethod <- match.arg(correlationMethod, c("pearson", "spearman", "kendall"))

    switch(
        EXPR = correlationMethod,
        pearson = .computePearsonCorrelations(assayMatrix),
        spearman = .computeSpearmanCorrelations(assayMatrix),
        kendall = .computeKendallCorrelations(assayMatrix)
    )
}

.reshapeCorrelationMatrixToEdges <- function(
    correlationMatrix,
    pValueMatrix,
    featureAnnotations,
    assayName,
    groupName,
    correlationMethod,
    correlationScope = "withinOmic",
    sourceType = "correlation",
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    pAdjustMethod = "fdr"
) {
    if (!length(correlationMatrix)) {
        return(.emptyStandardEdgeTable())
    }

    upperIndex <- upper.tri(correlationMatrix, diag = FALSE)
    fromIdentifiers <- rownames(correlationMatrix)[row(correlationMatrix)[upperIndex]]
    toIdentifiers <- colnames(correlationMatrix)[col(correlationMatrix)[upperIndex]]
    correlationValues <- as.numeric(correlationMatrix[upperIndex])
    pValues <- as.numeric(pValueMatrix[upperIndex])
    adjustedPValues <- stats::p.adjust(pValues, method = pAdjustMethod)

    edgeTable <- data.frame(
        fromFeatureIdentifier = fromIdentifiers,
        toFeatureIdentifier = toIdentifiers,
        correlationValue = correlationValues,
        pValue = pValues,
        adjustedPValue = adjustedPValues,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    keepIndex <- !is.na(edgeTable$correlationValue)
    if (!is.null(minimumAbsoluteCorrelation)) {
        keepIndex <- keepIndex & abs(edgeTable$correlationValue) >= minimumAbsoluteCorrelation
    }
    if (!is.null(adjustedPValueThreshold)) {
        keepIndex <- keepIndex & edgeTable$adjustedPValue <= adjustedPValueThreshold
    }
    edgeTable <- edgeTable[keepIndex, , drop = FALSE]

    if (!nrow(edgeTable)) {
        return(.emptyStandardEdgeTable())
    }

    featureAnnotations <- as.data.frame(featureAnnotations, stringsAsFactors = FALSE)
    fromMatch <- match(edgeTable$fromFeatureIdentifier, featureAnnotations$featureIdentifier)
    toMatch <- match(edgeTable$toFeatureIdentifier, featureAnnotations$featureIdentifier)
    edgeTable$fromFeatureName <- featureAnnotations$featureName[fromMatch]
    edgeTable$toFeatureName <- featureAnnotations$featureName[toMatch]
    edgeTable$fromAssayName <- assayName
    edgeTable$toAssayName <- assayName
    edgeTable$edgeType <- "correlation"
    edgeTable$edgeDirection <- ifelse(edgeTable$correlationValue >= 0, "positive", "negative")
    edgeTable$sourceType <- sourceType
    edgeTable$correlationScope <- correlationScope
    edgeTable$correlationMethod <- correlationMethod
    edgeTable$knowledgeSource <- NA_character_
    edgeTable$groupName <- groupName
    edgeTable$comparisonName <- NA_character_
    edgeTable$group1CorrelationValue <- NA_real_
    edgeTable$group2CorrelationValue <- NA_real_
    edgeTable$group1PValue <- NA_real_
    edgeTable$group2PValue <- NA_real_
    edgeTable$group1AdjustedPValue <- NA_real_
    edgeTable$group2AdjustedPValue <- NA_real_
    edgeTable$zScoreDifference <- NA_real_
    edgeTable$edgeWeight <- edgeTable$correlationValue
    edgeTable$evidenceScore <- NA_real_
    edgeTable$isDirected <- FALSE

    .coerceStandardEdgeTable(edgeTable)
}

.subsetAssayForCorrelation <- function(
    analysisData,
    assayName,
    sampleIds,
    featureSubset = NULL
) {
    featureSubset <- .resolveFeatureSubset(
        analysisData = analysisData,
        assayName = assayName,
        featureSubset = featureSubset
    )
    assayMatrix <- .extractAssayMatrix(analysisData, assayName)
    sampleIds <- intersect(sampleIds, colnames(assayMatrix))

    if (length(sampleIds) < 3L) {
        stop(
            "At least three samples are required to compute correlations for assay `",
            assayName, "`.",
            call. = FALSE
        )
    }

    assayMatrix[featureSubset, sampleIds, drop = FALSE]
}

#' Create a Group-Specific Correlation Network
#'
#' Create a within-omic correlation network for one assay and one sample
#' group.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param assayName Name of the assay to analyse.
#' @param groupColumn Sample metadata column used to define the group.
#' @param groupLevel Single group label to analyse.
#' @param sampleIds Optional explicit sample identifiers.
#' @param featureSubset Optional feature identifiers to retain.
#' @param correlationMethod Correlation method. One of `"pearson"`,
#'   `"spearman"`, or `"kendall"`.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation
#'   retained in the final network.
#' @param adjustedPValueThreshold Maximum adjusted p-value retained in
#'   the final network.
#' @param pAdjustMethod Multiple-testing correction method.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional name used when storing the result.
#' @param storeResult Whether to store the result in
#'   `metadata(analysisData)$cornetto`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @export
createCorrelationNetwork <- function(
    analysisData,
    assayName,
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL,
    featureSubset = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    pAdjustMethod = "fdr",
    featureNameColumn = "featureName",
    resultName = NULL,
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationMethod <- match.arg(correlationMethod)
    sampleIds <- .resolveGroupSamples(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevel = groupLevel,
        sampleIds = sampleIds
    )

    if (!is.null(featureSubset) && !length(featureSubset)) {
        edgeTable <- .emptyStandardEdgeTable()
        if (!storeResult) {
            return(edgeTable)
        }

        if (is.null(groupLevel)) {
            groupLabel <- .makeResultName(assayName, "customSamples")
        } else {
            groupLabel <- groupLevel
        }
        if (is.null(resultName)) {
            resultName <- .makeResultName(assayName, groupLabel, correlationMethod)
        }

        return(.storeCorNettoResults(
            analysisData = analysisData,
            slotName = "correlationNetworks",
            resultObject = edgeTable,
            resultName = resultName
        ))
    }

    assayMatrix <- .subsetAssayForCorrelation(
        analysisData = analysisData,
        assayName = assayName,
        sampleIds = sampleIds,
        featureSubset = featureSubset
    )
    featureAnnotations <- .extractFeatureAnnotations(
        analysisData = analysisData,
        assayName = assayName,
        featureNameColumn = featureNameColumn
    )
    featureAnnotations <- featureAnnotations[
        featureAnnotations$featureIdentifier %in% rownames(assayMatrix),
        ,
        drop = FALSE
    ]

    correlationResult <- .dispatchCorrelationMethod(
        assayMatrix = assayMatrix,
        correlationMethod = correlationMethod
    )

    if (is.null(groupLevel)) {
        groupLabel <- .makeResultName(assayName, "customSamples")
    } else {
        groupLabel <- groupLevel
    }

    edgeTable <- .reshapeCorrelationMatrixToEdges(
        correlationMatrix = correlationResult$correlationMatrix,
        pValueMatrix = correlationResult$pValueMatrix,
        featureAnnotations = featureAnnotations,
        assayName = assayName,
        groupName = groupLabel,
        correlationMethod = correlationMethod,
        correlationScope = "withinOmic",
        sourceType = "correlation",
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        adjustedPValueThreshold = adjustedPValueThreshold,
        pAdjustMethod = pAdjustMethod
    )

    if (!storeResult) {
        return(edgeTable)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName(assayName, groupLabel, correlationMethod)
    }

    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "correlationNetworks",
        resultObject = edgeTable,
        resultName = resultName
    )
}

#' Create Correlation Networks for Multiple Assays or Groups
#'
#' Wrapper around `createCorrelationNetwork()` for repeated dense
#' within-omic network construction.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param assayNames Optional assay names. When `NULL`, all assays are
#'   processed.
#' @param groupColumn Sample metadata column defining the groups.
#' @param groupLevels Optional group labels. When `NULL`, all levels in
#'   `groupColumn` are used.
#' @param sampleIdGroups Optional named list of explicit sample
#'   identifiers. When supplied, `groupColumn` and `groupLevels` are not
#'   required.
#' @param featureSubset Optional feature subset applied to all assays or
#'   supplied as a named list by assay.
#' @param correlationMethod Correlation method.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation
#'   retained.
#' @param adjustedPValueThreshold Maximum adjusted p-value retained.
#' @param pAdjustMethod Multiple-testing correction method.
#' @param featureNameColumn Row-data column containing display names.
#' @param storeResult Whether to store each generated network.
#'
#' @return A named list of edge tables or an updated
#'   `MultiAssayExperiment`.
#' @export
createCorrelationNetworks <- function(
    analysisData,
    assayNames = NULL,
    groupColumn = NULL,
    groupLevels = NULL,
    sampleIdGroups = NULL,
    featureSubset = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    pAdjustMethod = "fdr",
    featureNameColumn = "featureName",
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationMethod <- match.arg(correlationMethod)

    if (is.null(assayNames)) {
        assayNames <- names(MultiAssayExperiment::experiments(analysisData))
    }

    if (is.null(sampleIdGroups)) {
        .assertScalarCharacter(groupColumn, "groupColumn")
        sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData), stringsAsFactors = FALSE)
        if (is.null(groupLevels)) {
            groupLevels <- unique(as.character(sampleData[[groupColumn]]))
        }
        sampleIdGroups <- base::setNames(as.list(groupLevels), groupLevels)
        sampleIdGroups <- lapply(
            names(sampleIdGroups),
            function(groupLevel) {
                .resolveGroupSamples(
                    analysisData = analysisData,
                    groupColumn = groupColumn,
                    groupLevel = groupLevel
                )
            }
        )
        names(sampleIdGroups) <- groupLevels
    }

    resultList <- list()
    updatedAnalysisData <- analysisData
    for (assayName in assayNames) {
        assayFeatureSubset <- featureSubset
        if (is.list(featureSubset) && !is.null(featureSubset[[assayName]])) {
            assayFeatureSubset <- featureSubset[[assayName]]
        }

        for (groupName in names(sampleIdGroups)) {
            edgeTable <- createCorrelationNetwork(
                analysisData = updatedAnalysisData,
                assayName = assayName,
                sampleIds = sampleIdGroups[[groupName]],
                featureSubset = assayFeatureSubset,
                correlationMethod = correlationMethod,
                minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
                adjustedPValueThreshold = adjustedPValueThreshold,
                pAdjustMethod = pAdjustMethod,
                featureNameColumn = featureNameColumn,
                resultName = .makeResultName(assayName, groupName, correlationMethod),
                storeResult = FALSE
            )
            resultList[[.makeResultName(assayName, groupName, correlationMethod)]] <- edgeTable

            if (storeResult) {
                updatedAnalysisData <- .storeCorNettoResults(
                    analysisData = updatedAnalysisData,
                    slotName = "correlationNetworks",
                    resultObject = edgeTable,
                    resultName = .makeResultName(assayName, groupName, correlationMethod)
                )
            }
        }
    }

    if (storeResult) {
        return(updatedAnalysisData)
    }

    resultList
}
