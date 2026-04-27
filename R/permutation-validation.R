.validatePermutationCount <- function(nPermutations) {
    if (!is.numeric(nPermutations) || length(nPermutations) != 1L ||
        is.na(nPermutations) || nPermutations < 1L) {
        stop("`nPermutations` must be a positive integer.", call. = FALSE)
    }

    as.integer(nPermutations)
}

.empiricalPValue <- function(observedValue, nullValues, alternative) {
    nullValues <- as.numeric(nullValues)
    nullValues <- nullValues[!is.na(nullValues)]
    if (!length(nullValues) || is.na(observedValue)) {
        return(NA_real_)
    }

    extremeCount <- switch(
        EXPR = alternative,
        two.sided = sum(abs(nullValues) >= abs(observedValue)),
        greater = sum(nullValues >= observedValue),
        less = sum(nullValues <= observedValue)
    )

    (1 + extremeCount) / (1 + length(nullValues))
}

.emptySparseValidationTable <- function() {
    S4Vectors::DataFrame(
        fromFeatureIdentifier = character(),
        toFeatureIdentifier = character(),
        fromFeatureName = character(),
        toFeatureName = character(),
        fromAssayName = character(),
        toAssayName = character(),
        correlationScope = character(),
        correlationMethod = character(),
        knowledgeSource = character(),
        groupName = character(),
        observedCorrelationValue = numeric(),
        observedPValue = numeric(),
        observedAdjustedPValue = numeric(),
        empiricalPValue = numeric(),
        empiricalAdjustedPValue = numeric(),
        nullMeanCorrelation = numeric(),
        nullSdCorrelation = numeric(),
        validPermutations = integer(),
        nPermutations = integer(),
        nullModel = character(),
        alternative = character(),
        evidenceScore = numeric(),
        check.names = FALSE
    )
}

.computePermutedSparseCorrelation <- function(
    analysisData,
    edgeRow,
    sampleIds,
    correlationMethod,
    nullModel,
    assayPermutations = NULL
) {
    fromMatrix <- .extractAssayMatrix(analysisData, edgeRow$fromAssayName)
    toMatrix <- .extractAssayMatrix(analysisData, edgeRow$toAssayName)

    sharedSamples <- Reduce(
        intersect,
        list(sampleIds, colnames(fromMatrix), colnames(toMatrix))
    )
    if (length(sharedSamples) < 3L) {
        return(NA_real_)
    }

    x <- as.numeric(fromMatrix[edgeRow$fromFeatureIdentifier, sharedSamples, drop = TRUE])
    if (identical(nullModel, "permuteEndpointSamples")) {
        y <- as.numeric(toMatrix[edgeRow$toFeatureIdentifier, sharedSamples, drop = TRUE])
        completeIndex <- stats::complete.cases(x, y)
        if (sum(completeIndex) < 3L) {
            return(NA_real_)
        }
        y <- sample(y[completeIndex])
        x <- x[completeIndex]
    } else {
        permutedSamples <- assayPermutations[[edgeRow$toAssayName]][sharedSamples]
        y <- as.numeric(toMatrix[edgeRow$toFeatureIdentifier, permutedSamples, drop = TRUE])
        completeIndex <- stats::complete.cases(x, y)
        if (sum(completeIndex) < 3L) {
            return(NA_real_)
        }
        x <- x[completeIndex]
        y <- y[completeIndex]
    }

    tryCatch(
        suppressWarnings(stats::cor(x = x, y = y, method = correlationMethod)),
        error = function(...) NA_real_
    )
}

.makeAssaySamplePermutations <- function(analysisData, assayNames, sampleIds) {
    permutationList <- lapply(
        assayNames,
        function(assayName) {
            assayMatrix <- .extractAssayMatrix(analysisData, assayName)
            availableSamples <- intersect(sampleIds, colnames(assayMatrix))
            stats::setNames(sample(availableSamples), availableSamples)
        }
    )
    names(permutationList) <- assayNames
    permutationList
}

#' Validate Sparse Correlation Edges by Permutation
#'
#' Compute empirical p-values for prior-guided sparse correlation edges by
#' repeatedly breaking sample alignment and recomputing correlations.
#'
#' @inheritParams createSparseMultiOmicCorrelations
#' @param nPermutations Number of permutations.
#' @param nullModel Null model used to break sample alignment.
#'   `"permuteEndpointSamples"` permutes the second endpoint separately
#'   for each edge. `"permuteAssaySampleLabels"` permutes sample labels
#'   once per target assay per permutation, preserving target-assay
#'   covariance.
#' @param alternative Empirical tail used for p-value calculation.
#' @param empiricalPValueThreshold Optional maximum adjusted empirical
#'   p-value retained.
#' @param seed Optional random seed.
#'
#' @return A validation `DataFrame` containing observed correlations,
#'   asymptotic p-values, empirical p-values, and permutation summaries,
#'   or an updated `MultiAssayExperiment` when `storeResult = TRUE`.
#' @export
validateSparseCorrelationEdges <- function(
    analysisData,
    knowledgeNetwork,
    assayNames = NULL,
    correlationScope = c("all", "withinOmic", "crossOmic"),
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL,
    correlationMethod = c("pearson", "spearman", "kendall"),
    nPermutations = 1000L,
    nullModel = c("permuteEndpointSamples", "permuteAssaySampleLabels"),
    alternative = c("two.sided", "greater", "less"),
    minimumAbsoluteCorrelation = NULL,
    empiricalPValueThreshold = NULL,
    pAdjustMethod = "fdr",
    featureNameColumn = "featureName",
    seed = NULL,
    resultName = NULL,
    storeResult = FALSE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationScope <- match.arg(correlationScope)
    correlationMethod <- match.arg(correlationMethod)
    nullModel <- match.arg(nullModel)
    alternative <- match.arg(alternative)
    nPermutations <- .validatePermutationCount(nPermutations)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    sampleIds <- .resolveGroupSamples(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevel = groupLevel,
        sampleIds = sampleIds
    )
    groupName <- .resolveSparseGroupName(groupColumn, groupLevel, sampleIds)
    candidateEdges <- .selectSparseCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = knowledgeNetwork,
        assayNames = assayNames,
        correlationScope = correlationScope
    )

    observedEdges <- .computeSparsePairCorrelations(
        analysisData = analysisData,
        candidateEdges = candidateEdges,
        sampleIds = sampleIds,
        correlationMethod = correlationMethod,
        minimumAbsoluteCorrelation = NULL,
        adjustedPValueThreshold = NULL,
        pAdjustMethod = pAdjustMethod,
        featureNameColumn = featureNameColumn,
        groupName = groupName
    )

    if (!nrow(observedEdges)) {
        validationTable <- .emptySparseValidationTable()
        if (!storeResult) {
            return(validationTable)
        }
        if (is.null(resultName)) {
            resultName <- .makeResultName("sparseEdgeValidation", correlationScope, groupName, correlationMethod)
        }
        return(.storeCorNettoResults(analysisData, "validationResults", validationTable, resultName))
    }

    observedData <- as.data.frame(observedEdges, stringsAsFactors = FALSE, check.names = FALSE)
    edgeCount <- nrow(observedData)
    nullValues <- matrix(NA_real_, nrow = edgeCount, ncol = nPermutations)
    uniqueTargetAssays <- unique(observedData$toAssayName)

    for (permutationIndex in seq_len(nPermutations)) {
        assayPermutations <- NULL
        if (identical(nullModel, "permuteAssaySampleLabels")) {
            assayPermutations <- .makeAssaySamplePermutations(
                analysisData = analysisData,
                assayNames = uniqueTargetAssays,
                sampleIds = sampleIds
            )
        }

        for (edgeIndex in seq_len(edgeCount)) {
            nullValues[edgeIndex, permutationIndex] <- .computePermutedSparseCorrelation(
                analysisData = analysisData,
                edgeRow = observedData[edgeIndex, , drop = FALSE],
                sampleIds = sampleIds,
                correlationMethod = correlationMethod,
                nullModel = nullModel,
                assayPermutations = assayPermutations
            )
        }
    }

    empiricalPValues <- vapply(
        seq_len(edgeCount),
        function(edgeIndex) {
            .empiricalPValue(
                observedValue = observedData$correlationValue[[edgeIndex]],
                nullValues = nullValues[edgeIndex, ],
                alternative = alternative
            )
        },
        numeric(1L)
    )
    validPermutationCounts <- rowSums(!is.na(nullValues))

    validationTable <- data.frame(
        fromFeatureIdentifier = observedData$fromFeatureIdentifier,
        toFeatureIdentifier = observedData$toFeatureIdentifier,
        fromFeatureName = observedData$fromFeatureName,
        toFeatureName = observedData$toFeatureName,
        fromAssayName = observedData$fromAssayName,
        toAssayName = observedData$toAssayName,
        correlationScope = observedData$correlationScope,
        correlationMethod = observedData$correlationMethod,
        knowledgeSource = observedData$knowledgeSource,
        groupName = observedData$groupName,
        observedCorrelationValue = observedData$correlationValue,
        observedPValue = observedData$pValue,
        observedAdjustedPValue = observedData$adjustedPValue,
        empiricalPValue = empiricalPValues,
        empiricalAdjustedPValue = .adjustPValuesWithMissing(empiricalPValues, pAdjustMethod),
        nullMeanCorrelation = rowMeans(nullValues, na.rm = TRUE),
        nullSdCorrelation = apply(nullValues, 1L, stats::sd, na.rm = TRUE),
        validPermutations = validPermutationCounts,
        nPermutations = nPermutations,
        nullModel = nullModel,
        alternative = alternative,
        evidenceScore = observedData$evidenceScore,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    validationTable$nullMeanCorrelation[validationTable$validPermutations == 0L] <- NA_real_
    validationTable$nullSdCorrelation[validationTable$validPermutations < 2L] <- NA_real_

    keepIndex <- rep(TRUE, nrow(validationTable))
    if (!is.null(minimumAbsoluteCorrelation)) {
        keepIndex <- keepIndex & abs(validationTable$observedCorrelationValue) >= minimumAbsoluteCorrelation
    }
    if (!is.null(empiricalPValueThreshold)) {
        keepIndex <- keepIndex &
            !is.na(validationTable$empiricalAdjustedPValue) &
            validationTable$empiricalAdjustedPValue <= empiricalPValueThreshold
    }
    validationTable <- S4Vectors::DataFrame(validationTable[keepIndex, , drop = FALSE], check.names = FALSE)

    if (!storeResult) {
        return(validationTable)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName("sparseEdgeValidation", correlationScope, groupName, correlationMethod)
    }

    .storeCorNettoResults(analysisData, "validationResults", validationTable, resultName)
}

.permuteGroupLabels <- function(sampleData, groupColumn, groupLevels, blockColumn = NULL) {
    sampleData <- as.data.frame(sampleData, stringsAsFactors = FALSE, check.names = FALSE)
    labels <- as.character(sampleData[[groupColumn]])
    eligibleSamples <- rownames(sampleData)[labels %in% groupLevels]
    permutedLabels <- labels

    if (is.null(blockColumn)) {
        permutedLabels[match(eligibleSamples, rownames(sampleData))] <- sample(labels[match(eligibleSamples, rownames(sampleData))])
        return(permutedLabels)
    }

    if (!blockColumn %in% names(sampleData)) {
        stop("`blockColumn` was not found in `colData`: ", blockColumn, call. = FALSE)
    }

    blockValues <- as.character(sampleData[[blockColumn]])
    for (blockValue in unique(blockValues[match(eligibleSamples, rownames(sampleData))])) {
        blockSamples <- eligibleSamples[blockValues[match(eligibleSamples, rownames(sampleData))] %in% blockValue]
        blockIndex <- match(blockSamples, rownames(sampleData))
        permutedLabels[blockIndex] <- sample(labels[blockIndex])
    }

    permutedLabels
}

.setTemporaryGroupColumn <- function(analysisData, groupColumn, groupLabels) {
    sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData), stringsAsFactors = FALSE)
    sampleData[[groupColumn]] <- groupLabels
    sampleData <- S4Vectors::DataFrame(sampleData, check.names = FALSE)
    rownames(sampleData) <- rownames(MultiAssayExperiment::colData(analysisData))
    MultiAssayExperiment::colData(analysisData) <- sampleData
    analysisData
}

.emptyPermutationScoreTable <- function() {
    S4Vectors::DataFrame(
        nodeKey = character(),
        permutation = integer(),
        rawRewiringScore = numeric(),
        rootMeanSquareRewiringScore = numeric(),
        degreeMatchedZScore = numeric(),
        degreeMatchedScaledScore = numeric(),
        check.names = FALSE
    )
}

.scoreSingleDifferentialRewiring <- function(
    analysisData,
    groupColumn,
    groupLevels,
    assayName,
    candidateEdgeTable,
    featureSubset,
    correlationMethod,
    minimumAbsoluteCorrelation,
    adjustedPValueThreshold,
    pAdjustMethod,
    featureNameColumn,
    differenceAdjustedPValueThreshold,
    edgeWeightMethod
) {
    differentialResults <- testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevels = groupLevels,
        assayName = assayName,
        candidateEdgeTable = candidateEdgeTable,
        featureSubset = featureSubset,
        correlationMethod = correlationMethod,
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        adjustedPValueThreshold = adjustedPValueThreshold,
        pAdjustMethod = pAdjustMethod,
        featureNameColumn = featureNameColumn,
        storeResult = FALSE
    )
    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = differenceAdjustedPValueThreshold,
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        edgeWeightMethod = edgeWeightMethod
    )
    rewiringTable <- calculateRewiringScores(
        differentialCorrelationNetwork = differentialNetwork,
        storeResult = FALSE
    )

    list(
        differentialCorrelationTable = differentialResults,
        differentialCorrelationNetwork = differentialNetwork,
        rewiringTable = rewiringTable
    )
}

#' Validate Differential Rewiring Scores by Permutation
#'
#' Permute group labels, rerun differential correlation testing, rebuild
#' the differential network, and calculate empirical node-level rewiring
#' p-values.
#'
#' @inheritParams testDifferentialCorrelation
#' @inheritParams createDifferentialCorrelationNetwork
#' @param nPermutations Number of group-label permutations.
#' @param blockColumn Optional sample metadata column. When supplied,
#'   group labels are permuted within each block.
#' @param seed Optional random seed.
#' @param keepPermutationScores Whether to return the node-by-permutation
#'   score table.
#' @param resultName Optional storage name.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A named list containing the observed differential-correlation
#'   table, observed differential network, rewiring validation table, and
#'   optionally the permutation score table. When `storeResult = TRUE`,
#'   the list is stored in `metadata(analysisData)$cornetto$validationResults`
#'   and the updated `MultiAssayExperiment` is returned.
#' @export
permuteDifferentialRewiring <- function(
    analysisData,
    groupColumn,
    groupLevels,
    assayName = NULL,
    candidateEdgeTable = NULL,
    featureSubset = NULL,
    correlationMethod = c("pearson", "spearman"),
    minimumAbsoluteCorrelation = 0.3,
    adjustedPValueThreshold = 0.05,
    pAdjustMethod = "fdr",
    featureNameColumn = "featureName",
    differenceAdjustedPValueThreshold = 0.05,
    edgeWeightMethod = c("scaledAbsoluteZScore", "absoluteZScore", "signedZScore"),
    nPermutations = 1000L,
    blockColumn = NULL,
    seed = NULL,
    keepPermutationScores = FALSE,
    resultName = NULL,
    storeResult = FALSE
) {
    analysisData <- validateAnalysisData(analysisData)
    correlationMethod <- match.arg(correlationMethod)
    edgeWeightMethod <- match.arg(edgeWeightMethod)
    nPermutations <- .validatePermutationCount(nPermutations)
    if (length(groupLevels) != 2L) {
        stop("`groupLevels` must contain exactly two group labels.", call. = FALSE)
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }

    observedResult <- .scoreSingleDifferentialRewiring(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevels = groupLevels,
        assayName = assayName,
        candidateEdgeTable = candidateEdgeTable,
        featureSubset = featureSubset,
        correlationMethod = correlationMethod,
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        adjustedPValueThreshold = adjustedPValueThreshold,
        pAdjustMethod = pAdjustMethod,
        featureNameColumn = featureNameColumn,
        differenceAdjustedPValueThreshold = differenceAdjustedPValueThreshold,
        edgeWeightMethod = edgeWeightMethod
    )

    observedRewiring <- as.data.frame(
        observedResult$rewiringTable,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData), stringsAsFactors = FALSE)
    temporaryGroupColumn <- ".cornettoPermutationGroup"
    while (temporaryGroupColumn %in% names(sampleData)) {
        temporaryGroupColumn <- paste0(temporaryGroupColumn, "_")
    }

    permutationScoreList <- vector("list", nPermutations)
    for (permutationIndex in seq_len(nPermutations)) {
        permutedLabels <- .permuteGroupLabels(
            sampleData = sampleData,
            groupColumn = groupColumn,
            groupLevels = groupLevels,
            blockColumn = blockColumn
        )
        permutedAnalysisData <- .setTemporaryGroupColumn(
            analysisData = analysisData,
            groupColumn = temporaryGroupColumn,
            groupLabels = permutedLabels
        )

        permutationResult <- .scoreSingleDifferentialRewiring(
            analysisData = permutedAnalysisData,
            groupColumn = temporaryGroupColumn,
            groupLevels = groupLevels,
            assayName = assayName,
            candidateEdgeTable = candidateEdgeTable,
            featureSubset = featureSubset,
            correlationMethod = correlationMethod,
            minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
            adjustedPValueThreshold = adjustedPValueThreshold,
            pAdjustMethod = pAdjustMethod,
            featureNameColumn = featureNameColumn,
            differenceAdjustedPValueThreshold = differenceAdjustedPValueThreshold,
            edgeWeightMethod = edgeWeightMethod
        )

        if (!nrow(permutationResult$rewiringTable)) {
            next
        }

        permutationScores <- as.data.frame(
            permutationResult$rewiringTable,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        permutationScoreList[[permutationIndex]] <- data.frame(
            nodeKey = permutationScores$nodeKey,
            permutation = permutationIndex,
            rawRewiringScore = permutationScores$rawRewiringScore,
            rootMeanSquareRewiringScore = permutationScores$rootMeanSquareRewiringScore,
            degreeMatchedZScore = permutationScores$degreeMatchedZScore,
            degreeMatchedScaledScore = permutationScores$degreeMatchedScaledScore,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    permutationScores <- Filter(Negate(is.null), permutationScoreList)
    if (length(permutationScores)) {
        permutationScores <- do.call(rbind, permutationScores)
    } else {
        permutationScores <- as.data.frame(.emptyPermutationScoreTable(), stringsAsFactors = FALSE)
    }

    if (!nrow(observedRewiring)) {
        rewiringValidation <- observedResult$rewiringTable
    } else {
        empiricalPValues <- numeric(nrow(observedRewiring))
        nullMeans <- numeric(nrow(observedRewiring))
        nullSds <- numeric(nrow(observedRewiring))

        for (nodeIndex in seq_len(nrow(observedRewiring))) {
            nodeScores <- rep(0, nPermutations)
            nodePermutationScores <- permutationScores[
                permutationScores$nodeKey %in% observedRewiring$nodeKey[[nodeIndex]],
                ,
                drop = FALSE
            ]
            if (nrow(nodePermutationScores)) {
                nodeScores[nodePermutationScores$permutation] <- nodePermutationScores$rawRewiringScore
            }

            empiricalPValues[[nodeIndex]] <- .empiricalPValue(
                observedValue = observedRewiring$rawRewiringScore[[nodeIndex]],
                nullValues = nodeScores,
                alternative = "greater"
            )
            nullMeans[[nodeIndex]] <- mean(nodeScores)
            nullSds[[nodeIndex]] <- stats::sd(nodeScores)
        }

        observedRewiring$empiricalPValue <- empiricalPValues
        observedRewiring$empiricalAdjustedPValue <- .adjustPValuesWithMissing(empiricalPValues, pAdjustMethod)
        observedRewiring$nullMeanRawRewiringScore <- nullMeans
        observedRewiring$nullSdRawRewiringScore <- nullSds
        observedRewiring$nPermutations <- nPermutations
        observedRewiring$blockColumn <- if (is.null(blockColumn)) NA_character_ else blockColumn
        rewiringValidation <- S4Vectors::DataFrame(observedRewiring, check.names = FALSE)
    }

    result <- list(
        differentialCorrelationTable = observedResult$differentialCorrelationTable,
        differentialCorrelationNetwork = observedResult$differentialCorrelationNetwork,
        rewiringTable = rewiringValidation,
        permutationScores = if (isTRUE(keepPermutationScores)) {
            S4Vectors::DataFrame(permutationScores, check.names = FALSE)
        } else {
            NULL
        }
    )

    if (!storeResult) {
        return(result)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName(groupLevels[[1L]], "vs", groupLevels[[2L]], "rewiringPermutation", separator = "_")
    }

    .storeCorNettoResults(analysisData, "validationResults", result, resultName)
}
