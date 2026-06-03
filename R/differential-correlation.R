.boundCorrelation <- function(correlationValues, epsilon = 1e-06) {
    pmin(pmax(correlationValues, -1 + epsilon), 1 - epsilon)
}

.adjustDifferentialPValues <- function(pValues, pAdjustMethod) {
    adjustedValues <- rep(NA_real_, length(pValues))
    nonMissingIndex <- !is.na(pValues)
    if (!any(nonMissingIndex)) {
        return(adjustedValues)
    }

    observedPValues <- as.numeric(pValues[nonMissingIndex])
    if (identical(pAdjustMethod, "qvalue")) {
        if (!requireNamespace("qvalue", quietly = TRUE)) {
            stop(
                "The `qvalue` package is required when `pAdjustMethod = \"qvalue\"`.",
                call. = FALSE
            )
        }
        if (length(observedPValues) == 1L) {
            adjustedValues[nonMissingIndex] <- observedPValues
            return(adjustedValues)
        }

        qvalueResult <- tryCatch(
            qvalue::qvalue(observedPValues),
            error = function(...) NULL
        )
        if (is.null(qvalueResult)) {
            stop(
                "qvalue adjustment failed. Try `pAdjustMethod = \"fdr\"` or inspect the p-value distribution.",
                call. = FALSE
            )
        }

        adjustedValues[nonMissingIndex] <- qvalueResult$qvalues
        return(adjustedValues)
    }

    adjustedValues[nonMissingIndex] <- stats::p.adjust(observedPValues, method = pAdjustMethod)
    adjustedValues
}

.resolveDGCAColumn <- function(ddcorTable, groupLevel, suffix) {
    candidateColumns <- c(
        paste0(groupLevel, suffix),
        paste0(make.names(groupLevel), suffix)
    )
    matchedColumns <- intersect(candidateColumns, names(ddcorTable))
    if (!length(matchedColumns)) {
        stop(
            "DGCA output is missing the expected column for group `", groupLevel,
            "` and suffix `", suffix, "`.",
            call. = FALSE
        )
    }

    matchedColumns[[1L]]
}

.buildDifferentialCorrelationDesign <- function(group1SampleIds, group2SampleIds, groupLevels) {
    sampleGroups <- factor(
        c(
            rep(groupLevels[[1L]], length(group1SampleIds)),
            rep(groupLevels[[2L]], length(group2SampleIds))
        ),
        levels = groupLevels
    )
    designMatrix <- stats::model.matrix(~ 0 + sampleGroups)
    colnames(designMatrix) <- groupLevels
    designMatrix
}

.computeDifferentialCorrelationDGCA <- function(
    analysisData,
    assayName,
    sharedFeatures,
    group1SampleIds,
    group2SampleIds,
    groupLevels,
    correlationMethod,
    featureNameColumn
) {
    assayMatrix <- .extractAssayMatrix(analysisData, assayName)
    group1SampleIds <- intersect(group1SampleIds, colnames(assayMatrix))
    group2SampleIds <- intersect(group2SampleIds, colnames(assayMatrix))

    if (length(group1SampleIds) < 4L || length(group2SampleIds) < 4L) {
        stop(
            "Differential correlation requires at least four samples from each group ",
            "to be present in assay `", assayName, "`.",
            call. = FALSE
        )
    }

    assayMatrix <- assayMatrix[sharedFeatures, c(group1SampleIds, group2SampleIds), drop = FALSE]
    if (nrow(assayMatrix) < 2L) {
        return(.emptyStandardEdgeTable())
    }

    featureAnnotations <- .asPlainDataFrame(
        .extractFeatureAnnotations(
            analysisData = analysisData,
            assayName = assayName,
            featureNameColumn = featureNameColumn
        )
    )
    featureAnnotations <- featureAnnotations[
        featureAnnotations$featureIdentifier %in% rownames(assayMatrix),
        ,
        drop = FALSE
    ]

    designMatrix <- .buildDifferentialCorrelationDesign(
        group1SampleIds = group1SampleIds,
        group2SampleIds = group2SampleIds,
        groupLevels = groupLevels
    )
    dgcaResults <- DGCA::ddcorAll(
        inputMat = assayMatrix,
        design = designMatrix,
        compare = groupLevels,
        corrType = correlationMethod,
        nPairs = "all",
        sortBy = "zScoreDiff",
        adjust = "none",
        nPerms = 0,
        classify = FALSE,
        heatmapPlot = FALSE,
        verbose = FALSE
    )
    dgcaResults <- .asPlainDataFrame(dgcaResults)
    if (!nrow(dgcaResults)) {
        return(.emptyStandardEdgeTable())
    }

    group1CorrelationColumn <- .resolveDGCAColumn(dgcaResults, groupLevels[[1L]], "_cor")
    group2CorrelationColumn <- .resolveDGCAColumn(dgcaResults, groupLevels[[2L]], "_cor")
    group1PValueColumn <- .resolveDGCAColumn(dgcaResults, groupLevels[[1L]], "_pVal")
    group2PValueColumn <- .resolveDGCAColumn(dgcaResults, groupLevels[[2L]], "_pVal")

    rawResults <- data.frame(
        fromFeatureIdentifier = dgcaResults$Gene1,
        toFeatureIdentifier = dgcaResults$Gene2,
        fromFeatureName = featureAnnotations$featureName[
            match(dgcaResults$Gene1, featureAnnotations$featureIdentifier)
        ],
        toFeatureName = featureAnnotations$featureName[
            match(dgcaResults$Gene2, featureAnnotations$featureIdentifier)
        ],
        fromAssayName = assayName,
        toAssayName = assayName,
        edgeType = "differentialCorrelation",
        edgeDirection = NA_character_,
        sourceType = "differentialCorrelationTest",
        correlationScope = "withinOmic",
        correlationMethod = correlationMethod,
        knowledgeSource = NA_character_,
        groupName = NA_character_,
        comparisonName = NA_character_,
        correlationValue = NA_real_,
        group1CorrelationValue = as.numeric(dgcaResults[[group1CorrelationColumn]]),
        group2CorrelationValue = as.numeric(dgcaResults[[group2CorrelationColumn]]),
        pValue = as.numeric(dgcaResults$pValDiff),
        adjustedPValue = NA_real_,
        group1PValue = as.numeric(dgcaResults[[group1PValueColumn]]),
        group2PValue = as.numeric(dgcaResults[[group2PValueColumn]]),
        group1AdjustedPValue = NA_real_,
        group2AdjustedPValue = NA_real_,
        zScoreDifference = as.numeric(dgcaResults$zScoreDiff),
        sampleCount = NA_real_,
        group1SampleCount = length(group1SampleIds),
        group2SampleCount = length(group2SampleIds),
        edgeWeight = NA_real_,
        evidenceScore = NA_real_,
        isDirected = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    .coerceStandardEdgeTable(rawResults)
}

.computeDifferentialCorrelationFromCandidates <- function(
    analysisData,
    candidateEdgeTable,
    group1SampleIds,
    group2SampleIds,
    correlationMethod
) {
    candidateEdgeTable <- .coerceStandardEdgeTable(candidateEdgeTable)
    if (!nrow(candidateEdgeTable)) {
        return(.emptyStandardEdgeTable())
    }

    candidateData <- .asPlainDataFrame(candidateEdgeTable)
    availableAssays <- names(MultiAssayExperiment::experiments(analysisData))
    measuredIndex <- candidateData$fromAssayName %in% availableAssays &
        candidateData$toAssayName %in% availableAssays
    candidateData <- candidateData[measuredIndex, , drop = FALSE]
    if (!nrow(candidateData)) {
        return(.emptyStandardEdgeTable())
    }

    resultList <- vector("list", nrow(candidateData))

    for (edgeIndex in seq_len(nrow(candidateData))) {
        edgeRow <- candidateData[edgeIndex, , drop = FALSE]
        fromMatrix <- .extractAssayMatrix(analysisData, edgeRow$fromAssayName)
        toMatrix <- .extractAssayMatrix(analysisData, edgeRow$toAssayName)

        if (!edgeRow$fromFeatureIdentifier %in% rownames(fromMatrix)) {
            next
        }
        if (!edgeRow$toFeatureIdentifier %in% rownames(toMatrix)) {
            next
        }

        group1SharedSamples <- Reduce(
            intersect,
            list(group1SampleIds, colnames(fromMatrix), colnames(toMatrix))
        )
        group2SharedSamples <- Reduce(
            intersect,
            list(group2SampleIds, colnames(fromMatrix), colnames(toMatrix))
        )

        if (length(group1SharedSamples) < 4L || length(group2SharedSamples) < 4L) {
            next
        }

        group1X <- as.numeric(fromMatrix[edgeRow$fromFeatureIdentifier, group1SharedSamples, drop = TRUE])
        group1Y <- as.numeric(toMatrix[edgeRow$toFeatureIdentifier, group1SharedSamples, drop = TRUE])
        group2X <- as.numeric(fromMatrix[edgeRow$fromFeatureIdentifier, group2SharedSamples, drop = TRUE])
        group2Y <- as.numeric(toMatrix[edgeRow$toFeatureIdentifier, group2SharedSamples, drop = TRUE])

        group1Complete <- stats::complete.cases(group1X, group1Y)
        group2Complete <- stats::complete.cases(group2X, group2Y)
        if (sum(group1Complete) < 4L || sum(group2Complete) < 4L) {
            next
        }

        group1Arguments <- list(
            x = group1X[group1Complete],
            y = group1Y[group1Complete],
            method = correlationMethod
        )
        group2Arguments <- list(
            x = group2X[group2Complete],
            y = group2Y[group2Complete],
            method = correlationMethod
        )
        if (!identical(correlationMethod, "pearson")) {
            group1Arguments$exact <- FALSE
            group2Arguments$exact <- FALSE
        }

        group1Test <- tryCatch(
            suppressWarnings(do.call(stats::cor.test, group1Arguments)),
            error = function(...) NULL
        )
        group2Test <- tryCatch(
            suppressWarnings(do.call(stats::cor.test, group2Arguments)),
            error = function(...) NULL
        )

        if (is.null(group1Test) || is.null(group2Test)) {
            next
        }

        group1SampleCount <- sum(group1Complete)
        group2SampleCount <- sum(group2Complete)
        group1Correlation <- unname(group1Test$estimate[[1L]])
        group2Correlation <- unname(group2Test$estimate[[1L]])
        zScoreDifference <- DGCA::dCorrs(
            rho1 = .boundCorrelation(group1Correlation),
            n1 = group1SampleCount,
            rho2 = .boundCorrelation(group2Correlation),
            n2 = group2SampleCount,
            corrType = correlationMethod
        )
        pValueDifference <- 2 * stats::pnorm(abs(zScoreDifference), lower.tail = FALSE)

        resultList[[edgeIndex]] <- data.frame(
            fromFeatureIdentifier = edgeRow$fromFeatureIdentifier,
            toFeatureIdentifier = edgeRow$toFeatureIdentifier,
            fromFeatureName = edgeRow$fromFeatureName,
            toFeatureName = edgeRow$toFeatureName,
            fromAssayName = edgeRow$fromAssayName,
            toAssayName = edgeRow$toAssayName,
            edgeType = edgeRow$edgeType,
            edgeDirection = NA_character_,
            sourceType = "differentialCorrelationTest",
            correlationScope = if (!is.na(edgeRow$correlationScope) && nzchar(edgeRow$correlationScope)) {
                edgeRow$correlationScope
            } else if (identical(edgeRow$fromAssayName, edgeRow$toAssayName)) {
                "withinOmic"
            } else {
                "crossOmic"
            },
            correlationMethod = correlationMethod,
            knowledgeSource = edgeRow$knowledgeSource,
            groupName = NA_character_,
            comparisonName = NA_character_,
            correlationValue = NA_real_,
            group1CorrelationValue = group1Correlation,
            group2CorrelationValue = group2Correlation,
            pValue = pValueDifference,
            adjustedPValue = NA_real_,
            group1PValue = unname(group1Test$p.value),
            group2PValue = unname(group2Test$p.value),
            group1AdjustedPValue = NA_real_,
            group2AdjustedPValue = NA_real_,
            zScoreDifference = zScoreDifference,
            sampleCount = NA_real_,
            group1SampleCount = group1SampleCount,
            group2SampleCount = group2SampleCount,
            edgeWeight = NA_real_,
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

    .coerceStandardEdgeTable(do.call(rbind, resultList))
}

.finalizeDifferentialCorrelationResults <- function(
    edgeTable,
    groupLevels,
    minimumAbsoluteCorrelation,
    adjustedPValueThreshold,
    pAdjustMethod,
    correlationMethod
) {
    edgeTable <- .coerceStandardEdgeTable(edgeTable)
    if (!nrow(edgeTable)) {
        return(edgeTable)
    }

    edgeData <- .asPlainDataFrame(edgeTable)
    edgeData$comparisonName <- .makeResultName(groupLevels[[1L]], "vs", groupLevels[[2L]], separator = "_")
    edgeData$correlationMethod <- correlationMethod
    edgeData$sourceType <- "differentialCorrelationTest"
    edgeData$edgeType <- "differentialCorrelation"
    edgeData$edgeDirection <- NA_character_
    edgeData$isDirected <- FALSE

    group1Strong <- !is.na(edgeData$group1CorrelationValue)
    group2Strong <- !is.na(edgeData$group2CorrelationValue)
    if (!is.null(minimumAbsoluteCorrelation)) {
        group1Strong <- group1Strong &
            abs(edgeData$group1CorrelationValue) >= minimumAbsoluteCorrelation
        group2Strong <- group2Strong &
            abs(edgeData$group2CorrelationValue) >= minimumAbsoluteCorrelation
    }

    retainForAdjustment <- group1Strong | group2Strong
    edgeData <- edgeData[retainForAdjustment, , drop = FALSE]
    group1Strong <- group1Strong[retainForAdjustment]
    group2Strong <- group2Strong[retainForAdjustment]
    if (!nrow(edgeData)) {
        return(.emptyStandardEdgeTable())
    }

    edgeData$group1AdjustedPValue <- .adjustDifferentialPValues(edgeData$group1PValue, pAdjustMethod)
    edgeData$group2AdjustedPValue <- .adjustDifferentialPValues(edgeData$group2PValue, pAdjustMethod)
    edgeData$adjustedPValue <- .adjustDifferentialPValues(edgeData$pValue, pAdjustMethod)

    if (!is.null(adjustedPValueThreshold)) {
        group1Supported <- group1Strong &
            !is.na(edgeData$group1AdjustedPValue) &
            edgeData$group1AdjustedPValue <= adjustedPValueThreshold
        group2Supported <- group2Strong &
            !is.na(edgeData$group2AdjustedPValue) &
            edgeData$group2AdjustedPValue <= adjustedPValueThreshold
        keepIndex <- group1Supported | group2Supported
        edgeData <- edgeData[keepIndex, , drop = FALSE]
        if (!nrow(edgeData)) {
            return(.emptyStandardEdgeTable())
        }
    }

    .coerceStandardEdgeTable(edgeData)
}

.classifyCorrelationChange <- function(
    group1CorrelationValue,
    group2CorrelationValue,
    minimumAbsoluteCorrelation
) {
    presentInGroup1 <- abs(group1CorrelationValue) >= minimumAbsoluteCorrelation
    presentInGroup2 <- abs(group2CorrelationValue) >= minimumAbsoluteCorrelation

    if (presentInGroup1 && presentInGroup2 && sign(group1CorrelationValue) != sign(group2CorrelationValue)) {
        return("signFlip")
    }
    if (presentInGroup1 && !presentInGroup2) {
        return("gainInGroup1")
    }
    if (!presentInGroup1 && presentInGroup2) {
        return("lossInGroup1")
    }

    if (abs(group1CorrelationValue) > abs(group2CorrelationValue)) {
        if (group1CorrelationValue >= 0) {
            return("strengthenedPositive")
        }
        return("strengthenedNegative")
    }

    if (abs(group1CorrelationValue) < abs(group2CorrelationValue)) {
        if (group1CorrelationValue >= 0) {
            return("weakenedPositive")
        }
        return("weakenedNegative")
    }

    if (group1CorrelationValue >= 0) {
        return("positiveStable")
    }

    "negativeStable"
}

#' Test Differential Correlation Between Two Groups
#'
#' Compute differential correlation statistics for either all within-omic
#' feature pairs in one assay or a supplied set of candidate edges.
#' Within-omic all-pairs testing uses DGCA, while candidate-edge testing
#' evaluates only the supplied edges.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param groupColumn Sample metadata column used to define the two
#'   groups.
#' @param groupLevels A character vector of length two giving the group
#'   labels. The first group is treated as the reference for gain or loss
#'   labels.
#' @param assayName Assay name used when testing all within-omic pairs.
#' @param candidateEdgeTable Optional candidate edges to test. May be a
#'   standardized CorNetto edge table or a validated knowledge network.
#' @param featureSubset Optional feature subset used when `assayName` is
#'   supplied.
#' @param correlationMethod Correlation method. Differential correlation
#'   currently supports `"pearson"` and `"spearman"`.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation that
#'   must be present in at least one group before multiple-testing
#'   correction and edge retention.
#' @param adjustedPValueThreshold Maximum within-group adjusted p-value
#'   that must be reached in the same group as the correlation threshold
#'   before an edge is kept. Set `NULL` to skip within-group significance
#'   filtering after the correlation prefilter.
#' @param pAdjustMethod Multiple-testing correction method.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional name used when storing the result.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @details Edges are first filtered to pairs with sufficient absolute
#'   correlation in at least one group. Within-group and differential
#'   p-values are then adjusted on this retained edge universe. When
#'   `adjustedPValueThreshold` is not `NULL`, an edge is retained only if
#'   the same group is both sufficiently correlated and significant after
#'   multiple-testing correction. This rewiring-first behavior creates a
#'   stable scoring network without requiring every edge to pass a
#'   differential-correlation p-value threshold.
#' @examples
#' analysisData <- exampleAnalysisData()
#' differentialTable <- testDifferentialCorrelation(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
#'     storeResult = FALSE
#' )
#' differentialTable
#' @export
testDifferentialCorrelation <- function(
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
    resultName = NULL,
    storeResult = TRUE
) {
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    correlationMethod <- match.arg(correlationMethod)
    if (length(groupLevels) != 2L) {
        stop("`groupLevels` must contain exactly two group labels.", call. = FALSE)
    }

    group1SampleIds <- .resolveGroupSamples(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevel = groupLevels[[1L]]
    )
    group2SampleIds <- .resolveGroupSamples(
        analysisData = analysisData,
        groupColumn = groupColumn,
        groupLevel = groupLevels[[2L]]
    )
    if (length(group1SampleIds) < 10L || length(group2SampleIds) < 10L) {
        warning(
            "Differential correlation group sizes are below 10 samples for at least one group: ",
            groupLevels[[1L]], " = ", length(group1SampleIds), ", ",
            groupLevels[[2L]], " = ", length(group2SampleIds),
            ". Correlation estimates are unstable below 10 samples per group.",
            call. = FALSE
        )
    }

    if (is.null(candidateEdgeTable) && !is.null(featureSubset) && !length(featureSubset)) {
        rawResults <- .emptyStandardEdgeTable()
        if (!storeResult) {
            return(rawResults)
        }

        if (is.null(resultName)) {
            resultName <- .makeResultName(groupLevels[[1L]], "vs", groupLevels[[2L]], "raw", separator = "_")
        }

        return(.storeCorNettoResults(
            analysisData = analysisData,
            slotName = "differentialCorrelationResults",
            resultObject = rawResults,
            resultName = resultName
        ))
    }

    if (is.null(candidateEdgeTable)) {
        .assertScalarCharacter(assayName, "assayName")
        sharedFeatures <- .resolveFeatureSubset(
            analysisData = analysisData,
            assayName = assayName,
            featureSubset = featureSubset
        )
        rawResults <- .computeDifferentialCorrelationDGCA(
            analysisData = analysisData,
            assayName = assayName,
            sharedFeatures = sharedFeatures,
            group1SampleIds = group1SampleIds,
            group2SampleIds = group2SampleIds,
            groupLevels = groupLevels,
            correlationMethod = correlationMethod,
            featureNameColumn = featureNameColumn
        )
    } else {
        candidateEdgeTable <- .coerceStandardEdgeTable(candidateEdgeTable)
        rawResults <- .computeDifferentialCorrelationFromCandidates(
            analysisData = analysisData,
            candidateEdgeTable = candidateEdgeTable,
            group1SampleIds = group1SampleIds,
            group2SampleIds = group2SampleIds,
            correlationMethod = correlationMethod
        )
    }

    rawResults <- .finalizeDifferentialCorrelationResults(
        edgeTable = rawResults,
        groupLevels = groupLevels,
        minimumAbsoluteCorrelation = minimumAbsoluteCorrelation,
        adjustedPValueThreshold = adjustedPValueThreshold,
        pAdjustMethod = pAdjustMethod,
        correlationMethod = correlationMethod
    )

    if (!storeResult) {
        return(rawResults)
    }

    if (is.null(resultName)) {
        resultName <- .makeResultName(groupLevels[[1L]], "vs", groupLevels[[2L]], "raw", separator = "_")
    }

    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "differentialCorrelationResults",
        resultObject = rawResults,
        resultName = resultName
    )
}

#' Create a Weighted Differential Correlation Network
#'
#' Convert the output of `testDifferentialCorrelation()` into a weighted
#' edge table with rewiring-oriented labels.
#'
#' @param differentialCorrelationTable Output of
#'   `testDifferentialCorrelation()`.
#' @param differenceAdjustedPValueThreshold Optional maximum adjusted
#'   differential-correlation p-value retained in the network. The
#'   default `NULL` keeps all rewiring-scoring edges returned by
#'   `testDifferentialCorrelation()`.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation used to
#'   classify gains, losses, sign flips, and strengthening or weakening.
#' @param edgeWeightMethod Method used to generate the edge weight.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   weighted differential network.
#' @param resultName Name used when storing the weighted differential
#'   network.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' differentialTable <- testDifferentialCorrelation(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
#'     storeResult = FALSE
#' )
#' differentialNetwork <- createDifferentialCorrelationNetwork(
#'     differentialTable,
#'     minimumAbsoluteCorrelation = 0
#' )
#' differentialNetwork
#' @export
createDifferentialCorrelationNetwork <- function(
    differentialCorrelationTable,
    differenceAdjustedPValueThreshold = NULL,
    minimumAbsoluteCorrelation = 0.3,
    edgeWeightMethod = c("signedZScore", "absoluteZScore"),
    analysisData = NULL,
    resultName = "differentialCorrelationNetwork",
    storeResult = FALSE
) {
    edgeWeightMethod <- match.arg(edgeWeightMethod)
    differentialCorrelationTable <- .coerceStandardEdgeTable(differentialCorrelationTable)

    if (!nrow(differentialCorrelationTable)) {
        if (!storeResult) {
            return(differentialCorrelationTable)
        }
        .assertMultiAssayExperiment(analysisData)
        return(.storeCorNettoResults(
            analysisData,
            "differentialCorrelationNetworks",
            differentialCorrelationTable,
            resultName
        ))
    }

    edgeData <- .asPlainDataFrame(differentialCorrelationTable)
    if (!is.null(differenceAdjustedPValueThreshold)) {
        if (!is.numeric(differenceAdjustedPValueThreshold) ||
            length(differenceAdjustedPValueThreshold) != 1L ||
            is.na(differenceAdjustedPValueThreshold) ||
            differenceAdjustedPValueThreshold < 0) {
            stop(
                "`differenceAdjustedPValueThreshold` must be NULL or a non-negative numeric scalar.",
                call. = FALSE
            )
        }
        keepIndex <- !is.na(edgeData$adjustedPValue) &
            edgeData$adjustedPValue <= differenceAdjustedPValueThreshold
        edgeData <- edgeData[keepIndex, , drop = FALSE]
        if (!nrow(edgeData)) {
            edgeTable <- .emptyStandardEdgeTable()
            if (!storeResult) {
                return(edgeTable)
            }
            .assertMultiAssayExperiment(analysisData)
            return(.storeCorNettoResults(
                analysisData = analysisData,
                slotName = "differentialCorrelationNetworks",
                resultObject = edgeTable,
                resultName = resultName
            ))
        }
    }

    edgeData$edgeDirection <- vapply(
        seq_len(nrow(edgeData)),
        function(edgeIndex) {
            .classifyCorrelationChange(
                group1CorrelationValue = edgeData$group1CorrelationValue[[edgeIndex]],
                group2CorrelationValue = edgeData$group2CorrelationValue[[edgeIndex]],
                minimumAbsoluteCorrelation = minimumAbsoluteCorrelation
            )
        },
        character(1L)
    )

    if (identical(edgeWeightMethod, "absoluteZScore")) {
        edgeData$edgeWeight <- abs(edgeData$zScoreDifference)
    } else {
        edgeData$edgeWeight <- edgeData$zScoreDifference
    }

    edgeData$edgeType <- "differentialCorrelation"
    edgeData$sourceType <- "differentialCorrelation"
    edgeData$isDirected <- FALSE
    edgeTable <- .coerceStandardEdgeTable(edgeData)
    if (!storeResult) {
        return(edgeTable)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "differentialCorrelationNetworks",
        resultObject = edgeTable,
        resultName = resultName
    )
}

.computeDegreeMatchedScores <- function(rawRewiringScore, totalConnections) {
    degreeBins <- cut(
        x = log2(totalConnections + 1),
        breaks = seq(
            from = 0,
            to = max(log2(totalConnections + 1)) + 1,
            by = 1
        ),
        include.lowest = TRUE
    )
    zScores <- rawRewiringScore

    for (binName in unique(degreeBins)) {
        binIndex <- degreeBins == binName
        binValues <- rawRewiringScore[binIndex]
        if (length(unique(binValues)) == 1L) {
            zScores[binIndex] <- 0
        } else {
            zScores[binIndex] <- as.numeric(scale(binValues))
        }
    }

    list(
        degreeBin = degreeBins,
        degreeMatchedZScore = zScores
    )
}

.fillMissingRewiringNodeNames <- function(rewiringTable, nodeTable) {
    nodeNameMatch <- nodeTable$nodeName[match(rewiringTable$nodeKey, nodeTable$nodeKey)]
    useStoredName <- !is.na(nodeNameMatch) & nzchar(nodeNameMatch)
    rewiringTable$nodeName[useStoredName] <- nodeNameMatch[useStoredName]

    missingName <- is.na(rewiringTable$nodeName) | !nzchar(rewiringTable$nodeName)
    rewiringTable$nodeName[missingName] <- rewiringTable$nodeIdentifier[missingName]
    rewiringTable
}

#' Calculate Node Rewiring Scores
#'
#' Calculate raw, RMS-normalized, and degree-matched rewiring scores from
#' a weighted differential-correlation network.
#'
#' @param differentialCorrelationNetwork A standardized differential
#'   correlation edge table.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   rewiring scores.
#' @param resultName Name used when storing the rewiring scores.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A `DataFrame` of rewiring scores or an updated
#'   `MultiAssayExperiment`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' differentialTable <- testDifferentialCorrelation(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = c("Recovered", "PASC"),
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     pAdjustMethod = "fdr",
#'     storeResult = FALSE
#' )
#' differentialNetwork <- createDifferentialCorrelationNetwork(
#'     differentialTable,
#'     minimumAbsoluteCorrelation = 0
#' )
#' rewiringScores <- calculateRewiringScores(differentialNetwork)
#' rewiringScores
#' @export
calculateRewiringScores <- function(
    differentialCorrelationNetwork,
    analysisData = NULL,
    resultName = "rewiringScores",
    storeResult = FALSE
) {
    differentialCorrelationNetwork <- .coerceStandardEdgeTable(differentialCorrelationNetwork)
    if (!nrow(differentialCorrelationNetwork)) {
        rewiringTable <- S4Vectors::DataFrame(
            nodeKey = character(),
            nodeIdentifier = character(),
            nodeName = character(),
            assayName = character(),
            totalConnections = integer(),
            rawRewiringScore = numeric(),
            rootMeanSquareRewiringScore = numeric(),
            degreeBin = character(),
            degreeMatchedZScore = numeric(),
            check.names = FALSE
        )
        if (!storeResult) {
            return(rewiringTable)
        }

        .assertMultiAssayExperiment(analysisData)
        return(.storeCorNettoResults(analysisData, "rewiringResults", rewiringTable, resultName))
    }

    edgeData <- .asPlainDataFrame(differentialCorrelationNetwork)
    longTable <- rbind(
        data.frame(
            nodeKey = .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier),
            nodeIdentifier = edgeData$fromFeatureIdentifier,
            nodeName = edgeData$fromFeatureName,
            assayName = edgeData$fromAssayName,
            edgeWeight = edgeData$edgeWeight,
            stringsAsFactors = FALSE
        ),
        data.frame(
            nodeKey = .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier),
            nodeIdentifier = edgeData$toFeatureIdentifier,
            nodeName = edgeData$toFeatureName,
            assayName = edgeData$toAssayName,
            edgeWeight = edgeData$edgeWeight,
            stringsAsFactors = FALSE
        )
    )

    nodeKeys <- unique(longTable$nodeKey)
    rewiringTable <- lapply(
        nodeKeys,
        function(nodeKey) {
            nodeWeights <- longTable$edgeWeight[longTable$nodeKey %in% nodeKey]
            totalConnections <- length(nodeWeights)
            rawRewiringScore <- sqrt(sum(nodeWeights ^ 2, na.rm = TRUE))
            rootMeanSquareRewiringScore <- rawRewiringScore / sqrt(totalConnections)

            data.frame(
                nodeKey = nodeKey,
                nodeIdentifier = longTable$nodeIdentifier[match(nodeKey, longTable$nodeKey)],
                nodeName = longTable$nodeName[match(nodeKey, longTable$nodeKey)],
                assayName = longTable$assayName[match(nodeKey, longTable$nodeKey)],
                totalConnections = totalConnections,
                rawRewiringScore = rawRewiringScore,
                rootMeanSquareRewiringScore = rootMeanSquareRewiringScore,
                stringsAsFactors = FALSE
            )
        }
    )

    rewiringTable <- do.call(rbind, rewiringTable)
    degreeMatched <- .computeDegreeMatchedScores(
        rawRewiringScore = rewiringTable$rawRewiringScore,
        totalConnections = rewiringTable$totalConnections
    )
    storedNodeTable <- .asPlainDataFrame(
        .getStoredNodeTable(differentialCorrelationNetwork, fallbackToEdges = TRUE)
    )
    rewiringTable <- .fillMissingRewiringNodeNames(rewiringTable, storedNodeTable)
    rewiringTable$degreeBin <- as.character(degreeMatched$degreeBin)
    rewiringTable$degreeMatchedZScore <- as.numeric(degreeMatched$degreeMatchedZScore)
    rewiringTable <- S4Vectors::DataFrame(rewiringTable, check.names = FALSE)

    if (!storeResult) {
        return(rewiringTable)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "rewiringResults",
        resultObject = rewiringTable,
        resultName = resultName
    )
}
