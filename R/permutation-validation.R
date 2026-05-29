.validatePermutationCount <- function(nPermutations) {
    if (!is.numeric(nPermutations) || length(nPermutations) != 1L ||
        is.na(nPermutations) || nPermutations < 1L ||
        nPermutations != as.integer(nPermutations)) {
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

.permuteGroupLabels <- function(sampleData, groupColumn, groupLevels, blockColumn = NULL) {
    sampleData <- .asPlainDataFrame(sampleData)
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
    sampleData <- .asPlainDataFrame(MultiAssayExperiment::colData(analysisData))
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

#' Permute Rewiring Scores
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
#' @param scoreColumn Rewiring-score column used to compute empirical
#'   node-level p-values.
#' @param resultName Optional storage name.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A named list containing the observed differential-correlation
#'   table, observed differential network, rewiring validation table, and
#'   optionally the permutation score table. When `storeResult = TRUE`,
#'   the list is stored in `metadata(analysisData)$cornetto$validationResults`
#'   and the updated `MultiAssayExperiment` is returned.
#' @examples
#' analysisData <- exampleAnalysisData()
#' validation <- permuteRewiringScores(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
#'     differenceAdjustedPValueThreshold = 1,
#'     nPermutations = 2,
#'     seed = 1
#' )
#' names(validation)
#' @export
permuteRewiringScores <- function(
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
    edgeWeightMethod = c("absoluteZScore", "signedZScore"),
    scoreColumn = "rawRewiringScore",
    nPermutations = 1000L,
    blockColumn = NULL,
    seed = NULL,
    keepPermutationScores = FALSE,
    resultName = NULL,
    storeResult = FALSE
) {
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    correlationMethod <- match.arg(correlationMethod)
    edgeWeightMethod <- match.arg(edgeWeightMethod)
    .assertScalarCharacter(scoreColumn, "scoreColumn")
    nPermutations <- .validatePermutationCount(nPermutations)
    if (length(groupLevels) != 2L) {
        stop("`groupLevels` must contain exactly two group labels.", call. = FALSE)
    }

    if (!is.null(seed)) {
        withr::local_seed(seed)
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

    observedRewiring <- .asPlainDataFrame(observedResult$rewiringTable)
    sampleData <- .asPlainDataFrame(MultiAssayExperiment::colData(analysisData))
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

        permutationScores <- .asPlainDataFrame(permutationResult$rewiringTable)
        permutationScoreList[[permutationIndex]] <- data.frame(
            nodeKey = permutationScores$nodeKey,
            permutation = permutationIndex,
            rawRewiringScore = permutationScores$rawRewiringScore,
            rootMeanSquareRewiringScore = permutationScores$rootMeanSquareRewiringScore,
            degreeMatchedZScore = permutationScores$degreeMatchedZScore,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    permutationScores <- Filter(Negate(is.null), permutationScoreList)
    if (length(permutationScores)) {
        permutationScores <- do.call(rbind, permutationScores)
    } else {
        permutationScores <- .asPlainDataFrame(.emptyPermutationScoreTable())
    }

    if (!nrow(observedRewiring)) {
        rewiringValidation <- observedResult$rewiringTable
        rewiringValidation$empiricalPValue <- numeric()
        rewiringValidation$empiricalAdjustedPValue <- numeric()
        rewiringValidation$nullMeanScore <- numeric()
        rewiringValidation$nullSdScore <- numeric()
        rewiringValidation$scoreColumn <- character()
        rewiringValidation$nPermutations <- integer()
        rewiringValidation$blockColumn <- character()
    } else {
        if (!scoreColumn %in% names(observedRewiring)) {
            stop("`scoreColumn` was not found in the observed rewiring table.", call. = FALSE)
        }
        if (nrow(permutationScores) && !scoreColumn %in% names(permutationScores)) {
            stop("`scoreColumn` was not found in the permutation score table.", call. = FALSE)
        }
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
                nodeScores[nodePermutationScores$permutation] <- nodePermutationScores[[scoreColumn]]
            }

            empiricalPValues[[nodeIndex]] <- .empiricalPValue(
                observedValue = observedRewiring[[scoreColumn]][[nodeIndex]],
                nullValues = nodeScores,
                alternative = "greater"
            )
            nullMeans[[nodeIndex]] <- mean(nodeScores)
            nullSds[[nodeIndex]] <- stats::sd(nodeScores)
        }

        observedRewiring$empiricalPValue <- empiricalPValues
        observedRewiring$empiricalAdjustedPValue <- .adjustPValuesWithMissing(empiricalPValues, pAdjustMethod)
        observedRewiring$nullMeanScore <- nullMeans
        observedRewiring$nullSdScore <- nullSds
        observedRewiring$scoreColumn <- scoreColumn
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
