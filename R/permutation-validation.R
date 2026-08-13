.permutationTailProbability <- function(observedValue, nullValues, alternative) {
    nullValues <- as.numeric(nullValues)
    nullValues <- nullValues[!is.na(nullValues)]
    if (!length(nullValues) || is.na(observedValue)) {
        return(NA_real_)
    }

    extremeCount <- switch(
        EXPR = alternative,
        two.sided = sum(abs(nullValues) >= abs(observedValue)),
        greater = sum(nullValues >= observedValue),
        less = sum(nullValues <= observedValue),
        stop("Unknown `alternative`: ", alternative, call. = FALSE)
    )

    # The plus-one correction avoids a zero Monte Carlo estimate. When the
    # randomization design is valid, this is the Phipson-Smyth p-value.
    (1 + extremeCount) / (1 + length(nullValues))
}

.differentialEdgeKeys <- function(edgeTable) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    if (!nrow(edgeTable)) {
        return(character())
    }

    edgeData <- .asPlainDataFrame(edgeTable)
    fromKey <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    toKey <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)
    sort(paste(pmin(fromKey, toKey), pmax(fromKey, toKey), sep = "\r"))
}

.sameDifferentialEdgeSet <- function(first, second) {
    identical(.differentialEdgeKeys(first), .differentialEdgeKeys(second))
}

.candidateEdgesAlwaysEstimable <- function(
    analysisData,
    candidateEdgeTable,
    sampleData,
    groupColumn,
    groupLevels
) {
    labels <- as.character(sampleData[[groupColumn]])
    sampleIds <- rownames(sampleData)[labels %in% groupLevels]
    groupSizes <- vapply(
        groupLevels,
        function(groupLevel) sum(labels == groupLevel, na.rm = TRUE),
        integer(1L)
    )
    smallestGroup <- min(groupSizes)
    if (smallestGroup < 4L) {
        return(FALSE)
    }

    edgeData <- .asPlainDataFrame(candidateEdgeTable)
    usedAssays <- unique(c(edgeData$fromAssayName, edgeData$toAssayName))
    assayCache <- lapply(usedAssays, function(a) .extractAssayMatrix(analysisData, a))
    names(assayCache) <- usedAssays

    for (i in seq_len(nrow(edgeData))) {
        edge <- edgeData[i, , drop = FALSE]
        fromMatrix <- assayCache[[edge$fromAssayName]]
        toMatrix <- assayCache[[edge$toAssayName]]
        fromValues <- rep(NA_real_, length(sampleIds))
        toValues <- rep(NA_real_, length(sampleIds))

        fromIndex <- match(sampleIds, colnames(fromMatrix))
        toIndex <- match(sampleIds, colnames(toMatrix))
        fromPresent <- !is.na(fromIndex)
        toPresent <- !is.na(toIndex)
        fromValues[fromPresent] <- fromMatrix[
            edge$fromFeatureIdentifier,
            fromIndex[fromPresent],
            drop = TRUE
        ]
        toValues[toPresent] <- toMatrix[
            edge$toFeatureIdentifier,
            toIndex[toPresent],
            drop = TRUE
        ]

        complete <- is.finite(fromValues) & is.finite(toValues)
        missingCount <- sum(!complete)
        if (missingCount > smallestGroup - 4L) {
            return(FALSE)
        }

        largestFromTie <- max(tabulate(match(fromValues[complete], unique(fromValues[complete]))))
        largestToTie <- max(tabulate(match(toValues[complete], unique(toValues[complete]))))
        if (largestFromTie + missingCount >= smallestGroup ||
            largestToTie + missingCount >= smallestGroup) {
            return(FALSE)
        }
    }

    TRUE
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
    edgeWeightMethod,
    fixedEdgeUniverse = NULL
) {
    if (is.null(fixedEdgeUniverse)) {
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
    } else {
        # Permutations must score the edges the observed run scored. Refiltering
        # per permutation gives every permutation a different network, and the
        # rewiring score grows with degree, so the null would end up mixing
        # rewiring magnitude with network size.
        differentialResults <- testDifferentialCorrelation(
            analysisData = analysisData,
            groupColumn = groupColumn,
            groupLevels = groupLevels,
            candidateEdgeTable = fixedEdgeUniverse,
            correlationMethod = correlationMethod,
            minimumAbsoluteCorrelation = 0,
            adjustedPValueThreshold = NULL,
            pAdjustMethod = pAdjustMethod,
            featureNameColumn = featureNameColumn,
            storeResult = FALSE
        )
    }

    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = if (is.null(fixedEdgeUniverse)) {
            differenceAdjustedPValueThreshold
        } else {
            NULL
        },
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

.makePermutationJobs <- function(indices, permutedLabels, nWorkers) {
    nJobs <- min(length(indices), nWorkers)
    breaks <- floor(seq.int(0, length(indices), length.out = nJobs + 1L))

    lapply(seq_len(nJobs), function(i) {
        positions <- seq.int(breaks[[i]] + 1L, breaks[[i + 1L]])
        jobIndices <- indices[positions]
        list(
            indices = jobIndices,
            labels = permutedLabels[jobIndices]
        )
    })
}

.scorePermutationJob <- function(job, analysisData, settings) {
    scoreTables <- vector("list", length(job$indices))
    edgeSetComplete <- rep(TRUE, length(job$indices))

    for (i in seq_along(job$indices)) {
        permutationIndex <- job$indices[[i]]
        permutedAnalysisData <- .setTemporaryGroupColumn(
            analysisData = analysisData,
            groupColumn = settings$temporaryGroupColumn,
            groupLabels = job$labels[[i]]
        )

        # The observed run has already warned about group and assay sizes;
        # repeating that once per permutation is noise.
        permutationResult <- suppressWarnings(.scoreSingleDifferentialRewiring(
            analysisData = permutedAnalysisData,
            groupColumn = settings$temporaryGroupColumn,
            groupLevels = settings$groupLevels,
            assayName = settings$assayName,
            candidateEdgeTable = settings$candidateEdgeTable,
            featureSubset = settings$featureSubset,
            correlationMethod = settings$correlationMethod,
            minimumAbsoluteCorrelation = settings$minimumAbsoluteCorrelation,
            adjustedPValueThreshold = settings$adjustedPValueThreshold,
            pAdjustMethod = settings$pAdjustMethod,
            featureNameColumn = settings$featureNameColumn,
            differenceAdjustedPValueThreshold =
                settings$differenceAdjustedPValueThreshold,
            edgeWeightMethod = settings$edgeWeightMethod,
            fixedEdgeUniverse = settings$edgeUniverse
        ))

        if (settings$fixedCandidateDesign) {
            edgeSetComplete[[i]] <- .sameDifferentialEdgeSet(
                permutationResult$differentialCorrelationNetwork,
                settings$edgeUniverse
            )
        }

        if (!nrow(permutationResult$rewiringTable)) {
            next
        }

        permutationScores <- .asPlainDataFrame(permutationResult$rewiringTable)
        scoreTables[[i]] <- data.frame(
            nodeKey = permutationScores$nodeKey,
            permutation = permutationIndex,
            rawRewiringScore = permutationScores$rawRewiringScore,
            rootMeanSquareRewiringScore =
                permutationScores$rootMeanSquareRewiringScore,
            degreeMatchedZScore = permutationScores$degreeMatchedZScore,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    list(
        indices = job$indices,
        scoreTables = scoreTables,
        edgeSetComplete = edgeSetComplete
    )
}

#' Permute Rewiring Scores
#'
#' Permute group labels, rerun differential correlation testing, rebuild
#' the differential network, and compare observed node scores with their
#' permutation distributions. When both `assayName` and `candidateEdgeTable`
#' are `NULL`, every assay is tested densely within-omic; no cross-omic pairs
#' are generated.
#'
#' @inheritParams testDifferentialCorrelation
#' @inheritParams createDifferentialCorrelationNetwork
#' @param nPermutations Number of group-label permutations.
#' @param blockColumn Optional sample metadata column. When supplied,
#'   group labels are permuted within each block.
#' @param seed Optional random seed. Label allocations are generated before
#'   parallel work is dispatched, so the same seed gives the same allocations
#'   for serial and parallel backends.
#' @param keepPermutationScores Whether to return the node-by-permutation
#'   score table.
#' @param scoreColumn Rewiring-score column used to compute permutation tail
#'   probabilities.
#' @param resultName Optional storage name.
#' @param storeResult Whether to store the result in `analysisData`.
#' @param verbose Whether to report permutation progress.
#' @param progressEvery Number of completed permutations between progress
#'   messages. This also sets the size of each parallel batch.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   execution. The default [BiocParallel::SerialParam()] runs sequentially.
#'   On Windows, use [BiocParallel::SnowParam()] for parallel execution.
#'
#' @return A named list containing the observed differential-correlation
#'   table, observed differential network, rewiring validation table, and
#'   optionally the permutation score table. The rewiring table includes
#'   `permutationTailProbability`, `adjustedPermutationTailProbability`, and
#'   `inferenceStatus`; the same status is also a scalar element of the result
#'   list so it remains available when the table has no rows. When
#'   `storeResult = TRUE`,
#'   the list is stored in `metadata(analysisData)$cornetto$validationResults`
#'   and the updated `MultiAssayExperiment` is returned.
#' @details Group labels are permuted and the differential-correlation
#'   statistics are recomputed over a fixed edge universe. An
#'   edge can still be unscored in a permutation when a group has fewer than
#'   four pairwise-complete observations or a feature is constant, so realized
#'   node degree can fall in some permutations.
#'
#'   Dense all-pairs mode follows [testDifferentialCorrelation()]: with
#'   `assayName = NULL`, feature pairs are generated separately within every
#'   assay and the assay-specific edge sets are combined for scoring. It never
#'   creates cross-omic pairs. The resulting observed network defines the
#'   combined within-omic edge universe reused throughout permutation
#'   validation.
#'   Tail probabilities use the Phipson and Smyth (2010) estimator,
#'   `(1 + r) / (1 + B_i)`, where `B_i` is the number of permutations that
#'   produced a score for the node. `contributingPermutations` reports `B_i`.
#'   Nodes with few
#'   contributing permutations have a thin null and should be read as
#'   descriptive.
#'
#'   `inferenceStatus` is `"randomization p-value"` only when a supplied
#'   candidate universe is fixed before the observed analysis, no
#'   correlation or p-value filter is active, and missingness and ties cannot
#'   make a candidate edge unscorable under any permitted allocation. The code
#'   also checks that every edge was scored in the observed data and every
#'   sampled permutation. The selected node score must be finite in the
#'   observed data and all requested permutations. The user must still ensure
#'   that the candidate universe was prespecified and that labels are
#'   exchangeable globally, or within each level of `blockColumn`.
#'
#'   In all other settings, including all-pairs mode and the function's default
#'   filtering settings, `inferenceStatus` is `"conditional ranking"`. The two
#'   tail-probability columns then rank nodes against a label-dependent or
#'   incomplete reference distribution; they do not carry p-value or false
#'   discovery rate guarantees. Applying `pAdjustMethod` remains useful as a
#'   monotone ranking adjustment, but does not restore those guarantees.
#'
#'   Permutations are sampled independently and assignments can repeat. With
#'   `B = nPermutations` Monte Carlo draws, the smallest tail probability for a
#'   fully scored node is `1 / (B + 1)`. Exact enumeration of distinct assignments is not
#'   implemented. With 1,000 draws the minimum is about 0.001, which is often
#'   too coarse after adjustment across many nodes.
#'
#'   When `seed` is supplied, it is applied with [withr::local_seed()], which
#'   restores the calling RNG state on exit. With `seed = NULL`, the label
#'   draws use and advance the caller's current RNG stream.
#'   Label allocations are generated in permutation-index order before scoring
#'   begins. Parallel workers therefore perform deterministic calculations,
#'   and a fixed `seed` gives the same result regardless of worker count.
#'
#'   Parallel scoring is opt-in through `BPPARAM`; the default remains serial.
#'   When `verbose = TRUE`, the manager reports progress after each batch of up
#'   to `progressEvery` permutations. A stopped backend is started for the call
#'   and stopped on exit. A backend that was already running is left running.
#'   Each socket worker receives its own copy of the analysis inputs, so memory
#'   use should be considered when choosing the worker count. Interrupted runs
#'   cannot currently be resumed from a partial batch.
#' @references
#' Phipson B, Smyth GK (2010). Permutation p-values should never be zero:
#' calculating exact p-values when permutations are randomly drawn.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 9, 39.
#' \doi{10.2202/1544-6115.1585}
#' @examples
#' analysisData <- exampleAnalysisData()
#' validation <- permuteRewiringScores(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
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
    differenceAdjustedPValueThreshold = NULL,
    edgeWeightMethod = c("signedZScore", "absoluteZScore"),
    scoreColumn = "rawRewiringScore",
    nPermutations = 1000L,
    blockColumn = NULL,
    seed = NULL,
    keepPermutationScores = FALSE,
    resultName = NULL,
    storeResult = FALSE,
    verbose = FALSE,
    progressEvery = 300L,
    BPPARAM = BiocParallel::SerialParam()
) {
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    correlationMethod <- match.arg(correlationMethod)
    edgeWeightMethod <- match.arg(edgeWeightMethod)
    .assertScalarCharacter(scoreColumn, "scoreColumn")
    nPermutations <- .assertWholeNumber(nPermutations, "nPermutations", minimum = 1L)
    groupLevels <- .assertTwoGroupLevels(groupLevels)
    minimumAbsoluteCorrelation <- .validateUnitInterval(
        minimumAbsoluteCorrelation,
        "minimumAbsoluteCorrelation"
    )
    adjustedPValueThreshold <- .validateUnitInterval(
        adjustedPValueThreshold,
        "adjustedPValueThreshold",
        allowNull = TRUE
    )
    differenceAdjustedPValueThreshold <- .validateUnitInterval(
        differenceAdjustedPValueThreshold,
        "differenceAdjustedPValueThreshold",
        allowNull = TRUE
    )
    pAdjustMethod <- .validatePAdjustMethod(pAdjustMethod)
    .assertScalarCharacter(featureNameColumn, "featureNameColumn")
    .assertScalarLogical(keepPermutationScores, "keepPermutationScores")
    .assertScalarLogical(storeResult, "storeResult")
    .assertScalarLogical(verbose, "verbose")
    progressEvery <- .assertWholeNumber(
        progressEvery,
        "progressEvery",
        minimum = 1L
    )
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop(
            "`BPPARAM` must be a BiocParallelParam object.",
            call. = FALSE
        )
    }
    if (!is.null(blockColumn)) {
        .assertScalarCharacter(blockColumn, "blockColumn")
    }

    activeSelection <- (!is.null(minimumAbsoluteCorrelation) &&
        minimumAbsoluteCorrelation > 0) ||
        (!is.null(adjustedPValueThreshold) && adjustedPValueThreshold < 1) ||
        (!is.null(differenceAdjustedPValueThreshold) &&
            differenceAdjustedPValueThreshold < 1)
    if (activeSelection) {
        warning(
            "Observed-data edge filtering is active; permutation tail ",
            "probabilities are conditional rankings.",
            call. = FALSE
        )
    } else if (is.null(candidateEdgeTable)) {
        warning(
            "A prespecified `candidateEdgeTable` is required for randomization ",
            "inference; all-pairs permutation results are conditional rankings.",
            call. = FALSE
        )
    }

    if (!is.null(seed)) {
        withr::local_seed(seed)
    }

    fixedCandidateUniverse <- if (is.null(candidateEdgeTable)) {
        NULL
    } else {
        .prepareDifferentialCorrelationCandidates(analysisData, candidateEdgeTable)
    }
    observedCandidateEdges <- if (is.null(fixedCandidateUniverse)) {
        candidateEdgeTable
    } else {
        fixedCandidateUniverse
    }

    observedResult <- .scoreSingleDifferentialRewiring(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevels = groupLevels,
        assayName = assayName,
        candidateEdgeTable = observedCandidateEdges,
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

    fixedCandidateDesign <- !activeSelection &&
        !is.null(fixedCandidateUniverse) && nrow(fixedCandidateUniverse) > 0L
    edgeUniverse <- if (fixedCandidateDesign) {
        fixedCandidateUniverse
    } else {
        observedResult$differentialCorrelationNetwork
    }
    completeFixedUniverse <- fixedCandidateDesign &&
        .candidateEdgesAlwaysEstimable(
            analysisData = analysisData,
            candidateEdgeTable = fixedCandidateUniverse,
            sampleData = sampleData,
            groupColumn = groupColumn,
            groupLevels = groupLevels
        ) &&
        .sameDifferentialEdgeSet(
            observedResult$differentialCorrelationNetwork,
            edgeUniverse
        )

    permutedLabels <- lapply(seq_len(nPermutations), function(i) {
        .permuteGroupLabels(
            sampleData = sampleData,
            groupColumn = groupColumn,
            groupLevels = groupLevels,
            blockColumn = blockColumn
        )
    })

    permutationScoreList <- vector("list", nPermutations)
    edgeSetComplete <- rep(TRUE, nPermutations)
    nWorkers <- min(
        nPermutations,
        as.integer(BiocParallel::bpnworkers(BPPARAM))
    )
    settings <- list(
        temporaryGroupColumn = temporaryGroupColumn,
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
        edgeWeightMethod = edgeWeightMethod,
        edgeUniverse = edgeUniverse,
        fixedCandidateDesign = fixedCandidateDesign
    )

    backendWasRunning <- BiocParallel::bpisup(BPPARAM)
    if (!backendWasRunning) {
        BPPARAM <- BiocParallel::bpstart(BPPARAM)
        on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
    }

    if (verbose) {
        message(
            "Scoring ", nPermutations, " permutations with ", nWorkers,
            if (nWorkers == 1L) " worker." else " workers."
        )
    }

    batchStarts <- seq.int(1L, nPermutations, by = progressEvery)
    for (batchStart in batchStarts) {
        batchEnd <- min(batchStart + progressEvery - 1L, nPermutations)
        indices <- seq.int(batchStart, batchEnd)
        jobs <- .makePermutationJobs(
            indices = indices,
            permutedLabels = permutedLabels,
            nWorkers = nWorkers
        )
        jobResults <- BiocParallel::bplapply(
            jobs,
            .scorePermutationJob,
            analysisData = analysisData,
            settings = settings,
            BPPARAM = BPPARAM
        )

        for (jobResult in jobResults) {
            permutationScoreList[jobResult$indices] <- jobResult$scoreTables
            edgeSetComplete[jobResult$indices] <- jobResult$edgeSetComplete
        }

        if (verbose) {
            message(
                "Completed ", batchEnd, " of ", nPermutations,
                " permutations (",
                sprintf("%.1f", 100 * batchEnd / nPermutations),
                "%)."
            )
        }
    }
    completeFixedUniverse <- completeFixedUniverse && all(edgeSetComplete)

    permutationScores <- Filter(Negate(is.null), permutationScoreList)
    if (length(permutationScores)) {
        permutationScores <- do.call(rbind, permutationScores)
    } else {
        permutationScores <- data.frame(
            nodeKey = character(),
            permutation = integer(),
            rawRewiringScore = numeric(),
            rootMeanSquareRewiringScore = numeric(),
            degreeMatchedZScore = numeric(),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    scoreFullyEstimated <- FALSE
    if (!nrow(observedRewiring)) {
        rewiringValidation <- observedResult$rewiringTable
        rewiringValidation$permutationTailProbability <- numeric()
        rewiringValidation$adjustedPermutationTailProbability <- numeric()
        rewiringValidation$nullMeanScore <- numeric()
        rewiringValidation$nullSdScore <- numeric()
        rewiringValidation$contributingPermutations <- integer()
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
        tailProbabilities <- numeric(nrow(observedRewiring))
        nullMeans <- numeric(nrow(observedRewiring))
        nullSds <- numeric(nrow(observedRewiring))
        contributing <- integer(nrow(observedRewiring))

        for (i in seq_len(nrow(observedRewiring))) {
            # A node missing from a permutation is unscored, not a node scoring
            # zero. The rewiring score is a norm, so a zero can never count as
            # extreme, and padding with zeros biases the tail probability
            # downwards by roughly the node's dropout rate. Leave it NA and
            # report how many permutations contributed.
            nodeScores <- rep(NA_real_, nPermutations)
            nodePermutationScores <- permutationScores[
                permutationScores$nodeKey %in% observedRewiring$nodeKey[[i]],
                ,
                drop = FALSE
            ]
            if (nrow(nodePermutationScores)) {
                nodeScores[nodePermutationScores$permutation] <-
                    nodePermutationScores[[scoreColumn]]
            }

            contributing[[i]] <- sum(!is.na(nodeScores))
            tailProbabilities[[i]] <- .permutationTailProbability(
                observedValue = observedRewiring[[scoreColumn]][[i]],
                nullValues = nodeScores,
                alternative = "greater"
            )
            if (contributing[[i]] == 0L) {
                nullMeans[[i]] <- NA_real_
                nullSds[[i]] <- NA_real_
            } else {
                nullMeans[[i]] <- mean(nodeScores, na.rm = TRUE)
                nullSds[[i]] <- stats::sd(nodeScores, na.rm = TRUE)
            }
        }

        observedRewiring$permutationTailProbability <- tailProbabilities
        observedRewiring$adjustedPermutationTailProbability <-
            .adjustPValuesWithMissing(tailProbabilities, pAdjustMethod)
        observedRewiring$nullMeanScore <- nullMeans
        observedRewiring$nullSdScore <- nullSds
        observedRewiring$contributingPermutations <- contributing
        observedRewiring$scoreColumn <- scoreColumn
        observedRewiring$nPermutations <- nPermutations
        observedRewiring$blockColumn <- if (is.null(blockColumn)) NA_character_ else blockColumn
        scoreFullyEstimated <- all(is.finite(observedRewiring[[scoreColumn]])) &&
            all(contributing == nPermutations)
        completeFixedUniverse <- completeFixedUniverse && scoreFullyEstimated
        rewiringValidation <- S4Vectors::DataFrame(observedRewiring, check.names = FALSE)
    }

    inferenceStatus <- if (completeFixedUniverse) {
        "randomization p-value"
    } else {
        "conditional ranking"
    }
    rewiringValidation$inferenceStatus <- rep(inferenceStatus, nrow(rewiringValidation))

    if (fixedCandidateDesign && !completeFixedUniverse) {
        warning(
            "At least one prespecified edge or selected node score could not ",
            "be evaluated in the observed analysis and every permutation; ",
            "permutation tail probabilities are conditional rankings.",
            call. = FALSE
        )
    }

    result <- list(
        differentialCorrelationTable = observedResult$differentialCorrelationTable,
        differentialCorrelationNetwork = observedResult$differentialCorrelationNetwork,
        rewiringTable = rewiringValidation,
        inferenceStatus = inferenceStatus,
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
