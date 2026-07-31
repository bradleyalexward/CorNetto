# atanh() diverges at |r| = 1. The small bound keeps the Fisher transform
# finite for perfectly collinear pairs and is applied identically in dense and
# candidate-edge analyses.
.boundCorrelation <- function(correlationValues, epsilon = 1e-06) {
    pmin(pmax(correlationValues, -1 + epsilon), 1 - epsilon)
}

.fisherCorrelationDifference <- function(
    group1Correlation,
    group1SampleCount,
    group2Correlation,
    group2SampleCount,
    correlationMethod
) {
    varianceFactor <- if (identical(correlationMethod, "spearman")) 1.06 else 1
    numerator <- atanh(.boundCorrelation(group2Correlation)) -
        atanh(.boundCorrelation(group1Correlation))
    denominator <- sqrt(
        varianceFactor / (group1SampleCount - 3) +
            varianceFactor / (group2SampleCount - 3)
    )
    numerator / denominator
}

.computeDifferentialCorrelationAllPairs <- function(
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

    # The Fisher variance divides by (n - 3), so four per group is the arithmetic floor, not
    # a statistical recommendation. At 4 vs 4 the standard error of the z
    # difference is sqrt(2); .warnOnSmallGroups() is where that gets flagged.
    if (length(group1SampleIds) < 4L || length(group2SampleIds) < 4L) {
        stop(
            "Differential correlation requires at least four samples from each group ",
            "to be present in assay `", assayName, "`.",
            call. = FALSE
        )
    }
    .warnOnSmallGroups(groupLevels, length(group1SampleIds), length(group2SampleIds))

    assayMatrix <- assayMatrix[
        sharedFeatures,
        c(group1SampleIds, group2SampleIds),
        drop = FALSE
    ]
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

    group1Result <- .dispatchCorrelationMethod(
        assayMatrix = assayMatrix[, group1SampleIds, drop = FALSE],
        correlationMethod = correlationMethod
    )
    group2Result <- .dispatchCorrelationMethod(
        assayMatrix = assayMatrix[, group2SampleIds, drop = FALSE],
        correlationMethod = correlationMethod
    )

    upperIndex <- upper.tri(group1Result$correlationMatrix, diag = FALSE)
    group1Correlation <- as.numeric(group1Result$correlationMatrix[upperIndex])
    group2Correlation <- as.numeric(group2Result$correlationMatrix[upperIndex])
    group1SampleCount <- as.integer(group1Result$sampleCountMatrix[upperIndex])
    group2SampleCount <- as.integer(group2Result$sampleCountMatrix[upperIndex])
    testable <- !is.na(group1Correlation) & !is.na(group2Correlation) &
        !is.na(group1SampleCount) & !is.na(group2SampleCount) &
        group1SampleCount >= 4L & group2SampleCount >= 4L
    if (!any(testable)) {
        return(.emptyStandardEdgeTable())
    }

    fromIdentifiers <- rownames(group1Result$correlationMatrix)[
        row(group1Result$correlationMatrix)[upperIndex]
    ][testable]
    toIdentifiers <- colnames(group1Result$correlationMatrix)[
        col(group1Result$correlationMatrix)[upperIndex]
    ][testable]
    group1Correlation <- group1Correlation[testable]
    group2Correlation <- group2Correlation[testable]
    group1SampleCount <- group1SampleCount[testable]
    group2SampleCount <- group2SampleCount[testable]
    zScoreDifference <- .fisherCorrelationDifference(
        group1Correlation = group1Correlation,
        group1SampleCount = group1SampleCount,
        group2Correlation = group2Correlation,
        group2SampleCount = group2SampleCount,
        correlationMethod = correlationMethod
    )

    rawResults <- data.frame(
        fromFeatureIdentifier = fromIdentifiers,
        toFeatureIdentifier = toIdentifiers,
        fromFeatureName = featureAnnotations$featureName[
            match(fromIdentifiers, featureAnnotations$featureIdentifier)
        ],
        toFeatureName = featureAnnotations$featureName[
            match(toIdentifiers, featureAnnotations$featureIdentifier)
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
        group1CorrelationValue = group1Correlation,
        group2CorrelationValue = group2Correlation,
        pValue = 2 * stats::pnorm(abs(zScoreDifference), lower.tail = FALSE),
        adjustedPValue = NA_real_,
        group1PValue = as.numeric(group1Result$pValueMatrix[upperIndex])[testable],
        group2PValue = as.numeric(group2Result$pValueMatrix[upperIndex])[testable],
        group1AdjustedPValue = NA_real_,
        group2AdjustedPValue = NA_real_,
        zScoreDifference = zScoreDifference,
        sampleCount = NA_real_,
        group1SampleCount = group1SampleCount,
        group2SampleCount = group2SampleCount,
        edgeWeight = NA_real_,
        evidenceScore = NA_real_,
        isDirected = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    rawResults <- rawResults[
        order(abs(rawResults$zScoreDifference), decreasing = TRUE),
        ,
        drop = FALSE
    ]
    .coerceStandardEdgeTable(rawResults)
}

.prepareDifferentialCorrelationCandidates <- function(analysisData, candidateEdgeTable) {
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

    fromKey <- .nodeKey(
        candidateData$fromAssayName,
        candidateData$fromFeatureIdentifier
    )
    toKey <- .nodeKey(
        candidateData$toAssayName,
        candidateData$toFeatureIdentifier
    )
    selfLoop <- fromKey == toKey
    if (any(selfLoop)) {
        warning(
            "Dropping ", sum(selfLoop), " self-loop candidate edge(s); ",
            "a feature cannot be correlated with itself.",
            call. = FALSE
        )
        candidateData <- candidateData[!selfLoop, , drop = FALSE]
        fromKey <- fromKey[!selfLoop]
        toKey <- toKey[!selfLoop]
    }
    if (!nrow(candidateData)) {
        return(.emptyStandardEdgeTable())
    }
    candidateKey <- paste(pmin(fromKey, toKey), pmax(fromKey, toKey), sep = "\r")
    candidateData <- candidateData[!duplicated(candidateKey), , drop = FALSE]

    experimentList <- MultiAssayExperiment::experiments(analysisData)
    featureIds <- lapply(experimentList, rownames)
    endpointsAvailable <- vapply(
        seq_len(nrow(candidateData)),
        function(i) {
            candidateData$fromFeatureIdentifier[[i]] %in%
                featureIds[[candidateData$fromAssayName[[i]]]] &&
                candidateData$toFeatureIdentifier[[i]] %in%
                    featureIds[[candidateData$toAssayName[[i]]]]
        },
        logical(1L)
    )
    candidateData <- candidateData[endpointsAvailable, , drop = FALSE]
    if (!nrow(candidateData)) {
        return(.emptyStandardEdgeTable())
    }

    .coerceStandardEdgeTable(candidateData)
}

.computeDifferentialCorrelationFromCandidates <- function(
    analysisData,
    candidateEdgeTable,
    group1SampleIds,
    group2SampleIds,
    correlationMethod
) {
    candidateEdgeTable <- .prepareDifferentialCorrelationCandidates(
        analysisData,
        candidateEdgeTable
    )
    if (!nrow(candidateEdgeTable)) {
        return(.emptyStandardEdgeTable())
    }
    candidateData <- .asPlainDataFrame(candidateEdgeTable)

    # Extraction revalidates the whole matrix, so pull each assay once instead
    # of twice per candidate edge.
    usedAssays <- unique(c(candidateData$fromAssayName, candidateData$toAssayName))
    assayCache <- lapply(usedAssays, function(a) .extractAssayMatrix(analysisData, a))
    names(assayCache) <- usedAssays

    out <- vector("list", nrow(candidateData))

    for (i in seq_len(nrow(candidateData))) {
        edge <- candidateData[i, , drop = FALSE]
        fromMatrix <- assayCache[[edge$fromAssayName]]
        toMatrix <- assayCache[[edge$toAssayName]]

        if (!edge$fromFeatureIdentifier %in% rownames(fromMatrix)) {
            next
        }
        if (!edge$toFeatureIdentifier %in% rownames(toMatrix)) {
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

        # The Fisher variance divides by (n - 3), so four per group is the arithmetic floor.
        # It is not a statistical recommendation; see the warning at ten in
        # testDifferentialCorrelation().
        if (length(group1SharedSamples) < 4L || length(group2SharedSamples) < 4L) {
            next
        }

        group1Test <- .runCorrelationTest(
            x = as.numeric(fromMatrix[edge$fromFeatureIdentifier, group1SharedSamples, drop = TRUE]),
            y = as.numeric(toMatrix[edge$toFeatureIdentifier, group1SharedSamples, drop = TRUE]),
            correlationMethod = correlationMethod,
            minimumSampleCount = 4L
        )
        group2Test <- .runCorrelationTest(
            x = as.numeric(fromMatrix[edge$fromFeatureIdentifier, group2SharedSamples, drop = TRUE]),
            y = as.numeric(toMatrix[edge$toFeatureIdentifier, group2SharedSamples, drop = TRUE]),
            correlationMethod = correlationMethod,
            minimumSampleCount = 4L
        )
        if (is.null(group1Test) || is.null(group2Test)) {
            next
        }

        zScoreDifference <- .fisherCorrelationDifference(
            group1Correlation = group1Test$correlationValue,
            group1SampleCount = group1Test$sampleCount,
            group2Correlation = group2Test$correlationValue,
            group2SampleCount = group2Test$sampleCount,
            correlationMethod = correlationMethod
        )
        pValueDifference <- 2 * stats::pnorm(abs(zScoreDifference), lower.tail = FALSE)

        out[[i]] <- data.frame(
            fromFeatureIdentifier = edge$fromFeatureIdentifier,
            toFeatureIdentifier = edge$toFeatureIdentifier,
            fromFeatureName = edge$fromFeatureName,
            toFeatureName = edge$toFeatureName,
            fromAssayName = edge$fromAssayName,
            toAssayName = edge$toAssayName,
            edgeType = edge$edgeType,
            edgeDirection = NA_character_,
            sourceType = "differentialCorrelationTest",
            correlationScope = if (!is.na(edge$correlationScope) && nzchar(edge$correlationScope)) {
                edge$correlationScope
            } else if (identical(edge$fromAssayName, edge$toAssayName)) {
                "withinOmic"
            } else {
                "crossOmic"
            },
            correlationMethod = correlationMethod,
            knowledgeSource = edge$knowledgeSource,
            groupName = NA_character_,
            comparisonName = NA_character_,
            correlationValue = NA_real_,
            group1CorrelationValue = group1Test$correlationValue,
            group2CorrelationValue = group2Test$correlationValue,
            pValue = pValueDifference,
            adjustedPValue = NA_real_,
            group1PValue = group1Test$pValue,
            group2PValue = group2Test$pValue,
            group1AdjustedPValue = NA_real_,
            group2AdjustedPValue = NA_real_,
            zScoreDifference = zScoreDifference,
            sampleCount = NA_real_,
            group1SampleCount = group1Test$sampleCount,
            group2SampleCount = group2Test$sampleCount,
            edgeWeight = NA_real_,
            evidenceScore = edge$evidenceScore,
            isDirected = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    out <- Filter(Negate(is.null), out)
    if (!length(out)) {
        return(.emptyStandardEdgeTable())
    }

    .coerceStandardEdgeTable(do.call(rbind, out))
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

    # Adjustment must cover the full set of tested hypotheses. Filtering first
    # on |r| would select on a quantity tied directly to the p-value and make
    # the resulting BH values anti-conservative.
    edgeData$group1AdjustedPValue <- .adjustPValuesWithMissing(
        edgeData$group1PValue,
        pAdjustMethod
    )
    edgeData$group2AdjustedPValue <- .adjustPValuesWithMissing(
        edgeData$group2PValue,
        pAdjustMethod
    )
    edgeData$adjustedPValue <- .adjustPValuesWithMissing(
        edgeData$pValue,
        pAdjustMethod
    )

    retainForScoring <- group1Strong | group2Strong
    edgeData <- edgeData[retainForScoring, , drop = FALSE]
    group1Strong <- group1Strong[retainForScoring]
    group2Strong <- group2Strong[retainForScoring]
    if (!nrow(edgeData)) {
        return(.emptyStandardEdgeTable())
    }

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
    # A group with no variance in one feature yields no correlation at all,
    # which is not the same as a correlation of zero. Refuse to label it rather
    # than letting NA reach the `if` conditions below.
    if (is.na(group1CorrelationValue) || is.na(group2CorrelationValue)) {
        return(NA_character_)
    }

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
#' All-pairs testing evaluates every within-omic pair, while candidate-edge
#' testing evaluates only the supplied edges.
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
#'   Candidate self-loops are omitted with a warning because a feature's
#'   correlation with itself is not an informative differential edge.
#' @param featureSubset Optional feature subset used when `assayName` is
#'   supplied.
#' @param correlationMethod Correlation method. Differential correlation
#'   currently supports `"pearson"` and `"spearman"`.
#' @param minimumAbsoluteCorrelation Minimum absolute correlation that
#'   must be present in at least one group for edge retention. P-values are
#'   adjusted over the full testable family before this filter is applied.
#' @param adjustedPValueThreshold Maximum within-group adjusted p-value
#'   that must be reached in the same group as the correlation threshold
#'   before an edge is kept. Set `NULL` to skip within-group significance
#'   filtering after the correlation prefilter.
#' @param pAdjustMethod Multiple-testing correction method. Use one of
#'   `stats::p.adjust.methods` or `"qvalue"`; the latter requires the suggested
#'   `qvalue` package.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional name used when storing the result.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @details Within-group and differential p-values are adjusted over all
#'   testable pairs before any correlation-strength filter is applied. Edges
#'   are then filtered to pairs with sufficient absolute correlation in at
#'   least one group. When
#'   `adjustedPValueThreshold` is not `NULL`, an edge is retained only if
#'   the same group is both sufficiently correlated and significant after
#'   multiple-testing correction. Differential-correlation adjusted p-values
#'   are reported but are not used for retention unless requested later by
#'   [createDifferentialCorrelationNetwork()].
#'
#'   Candidate edges are treated as undirected hypotheses. Self-loops are
#'   removed with a warning, and duplicate or mirrored rows are tested once,
#'   using the first occurrence.
#'
#'   The statistic is the Fisher z difference between two independent
#'   correlations, evaluated against a standard normal:
#'
#'   \deqn{z = (atanh(r2) - atanh(r1)) / sqrt(1/(n1 - 3) + 1/(n2 - 3))}
#'
#'   For `correlationMethod = "spearman"` the variance term is inflated by 1.06
#'   following Fieller, Hartley and Pearson (1957).
#'
#'   Sample correlations of exactly -1 or 1 make the Fisher transform
#'   infinite. CorNetto bounds them to `-1 + 1e-6` or `1 - 1e-6` before the
#'   transform. This keeps the statistic finite but makes results for perfectly
#'   collinear pairs sensitive to that numerical convention.
#'
#'   Mind the direction. The subtraction is group 2 minus group 1, so a
#'   **positive** `zScoreDifference` means the pair is more strongly correlated
#'   in `groupLevels[2]`. This is the opposite orientation to the
#'   `edgeDirection` labels from [createDifferentialCorrelationNetwork()],
#'   which are named from group 1's point of view (`gainInGroup1` marks a pair
#'   correlated in group 1 and not in group 2, and therefore carries a negative
#'   z). The magnitude, the p-value, and every rewiring score are unaffected by
#'   the convention, since they square or take the absolute value of `z`; only
#'   `edgeWeight` under `edgeWeightMethod = "signedZScore"` exposes the sign.
#'
#'   This carries assumptions that omic data routinely violate, and the
#'   package does not test them for you:
#'
#'   - Pearson z tests assume each pair is approximately bivariate normal.
#'     Log-abundance proteomics and metabolomics are often heavy-tailed or
#'     left-censored, where the approximation is anticonservative. Use
#'     `correlationMethod = "spearman"` for heavy-tailed, monotone
#'     non-linear, or rank-robust analyses; Pearson remains the default
#'     because it is conventional for normalized log-scale omics.
#'   - The two groups must be independent. Repeated measurements on the same
#'     subject in both groups break this.
#'   - The variance term divides by `n - 3`, so four samples per group is the
#'     arithmetic floor and the function stops below it. That floor is not a
#'     recommendation: at four per group the standard error of the z
#'     difference is `sqrt(2)`. A warning is emitted below ten samples per
#'     group, and results there should be read as descriptive.
#'
#'   Correlations use pairwise complete observations, so each edge carries its
#'   own effective sample size in `group1SampleCount` and `group2SampleCount`.
#'   Pairs where a feature is constant within a group yield no correlation and
#'   are dropped.
#' @references
#' Fisher RA (1921). On the "probable error" of a coefficient of correlation
#' deduced from a small sample. \emph{Metron}, 1, 3-32.
#'
#' Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation
#' coefficients I. \emph{Biometrika}, 44, 470-481.
#' \doi{10.1093/biomet/44.3-4.470}
#'
#' Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate.
#' \emph{Journal of the Royal Statistical Society B}, 57, 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
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
    groupLevels <- .assertTwoGroupLevels(groupLevels)
    minimumAbsoluteCorrelation <- .validateUnitInterval(
        minimumAbsoluteCorrelation,
        "minimumAbsoluteCorrelation",
        allowNull = TRUE
    )
    adjustedPValueThreshold <- .validateUnitInterval(
        adjustedPValueThreshold,
        "adjustedPValueThreshold",
        allowNull = TRUE
    )
    pAdjustMethod <- .validatePAdjustMethod(pAdjustMethod)
    .assertScalarCharacter(featureNameColumn, "featureNameColumn")
    .assertScalarLogical(storeResult, "storeResult")

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
    allAssays <- is.null(assayName)
    if (allAssays) {
        assayNames <- names(MultiAssayExperiment::experiments(analysisData))
    } else {
        .assertScalarCharacter(assayName, "assayName")
        assayNames <- assayName
    }

    rawResultList <- vector("list", length(assayNames))
    names(rawResultList) <- assayNames

    for (currentAssay in assayNames) {
        assayFeatureSubset <- featureSubset

        if (allAssays && !is.null(featureSubset)) {
            assayFeatureSubset <- .resolveAssayFilterArgument(
                argumentValue = featureSubset,
                assayName = currentAssay,
                default = NULL
            )

            if (!is.null(assayFeatureSubset)) {
                assayFeatureSubset <- unique(as.character(assayFeatureSubset))

                if (anyNA(assayFeatureSubset) ||
                    any(!nzchar(assayFeatureSubset))) {
                    stop(
                        "`featureSubset` must be non-missing and non-empty.",
                        call. = FALSE
                    )
                }

                assayFeatureSubset <- intersect(
                    rownames(.extractAssayMatrix(
                        analysisData,
                        currentAssay
                    )),
                    assayFeatureSubset
                )
            }
        }

        sharedFeatures <- .resolveFeatureSubset(
            analysisData = analysisData,
            assayName = currentAssay,
            featureSubset = assayFeatureSubset
        )

        rawResultList[[currentAssay]] <- if (length(sharedFeatures) < 2L) {
            .emptyStandardEdgeTable()
        } else {
            .computeDifferentialCorrelationAllPairs(
                analysisData = analysisData,
                assayName = currentAssay,
                sharedFeatures = sharedFeatures,
                group1SampleIds = group1SampleIds,
                group2SampleIds = group2SampleIds,
                groupLevels = groupLevels,
                correlationMethod = correlationMethod,
                featureNameColumn = featureNameColumn
            )
        }
    }

    rawResults <- if (length(rawResultList) == 1L) {
        rawResultList[[1L]]
    } else {
        combineNetworks(
            networkList = rawResultList,
            includeReverseEdges = FALSE,
            removeDuplicates = FALSE
        )
    }
} else {
        .warnOnSmallGroups(groupLevels, length(group1SampleIds), length(group2SampleIds))
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
#'   classify gains, losses, sign flips, and strengthening or weakening. Note
#'   that this is independent of the prefilter applied by
#'   [testDifferentialCorrelation()]: if you change one, change the other, or
#'   the labels will disagree with which edges are present.
#' @param edgeWeightMethod Method used to generate the edge weight.
#'   `"signedZScore"` keeps the sign of `zScoreDifference`; `"absoluteZScore"`
#'   takes its magnitude. See Details for what the sign means.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   weighted differential network.
#' @param resultName Name used when storing the weighted differential
#'   network.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized edge `DataFrame`, or the updated
#'   `MultiAssayExperiment` when `storeResult = TRUE`.
#' @details `edgeWeightMethod` does not affect any rewiring score, because
#'   [calculateRewiringScores()] squares the weight. It matters only when you
#'   read or export `edgeWeight` directly.
#'
#'   Under `"signedZScore"` the weight inherits the Fisher statistic's
#'   group-2-minus-group-1 orientation: a **positive**
#'   weight means the pair is more strongly correlated in the second of the
#'   `groupLevels` passed to [testDifferentialCorrelation()]. The
#'   `edgeDirection` labels run the other way, being named from group 1's point
#'   of view, so `gainInGroup1` edges carry negative weights.
#'
#'   `edgeDirection` is `NA` when either group has no correlation to compare,
#'   which happens when a feature is constant within that group.
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
    differenceAdjustedPValueThreshold <- .validateUnitInterval(
        differenceAdjustedPValueThreshold,
        "differenceAdjustedPValueThreshold",
        allowNull = TRUE
    )
    minimumAbsoluteCorrelation <- .validateUnitInterval(
        minimumAbsoluteCorrelation,
        "minimumAbsoluteCorrelation"
    )
    .assertScalarLogical(storeResult, "storeResult")
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
        function(i) {
            .classifyCorrelationChange(
                group1CorrelationValue = edgeData$group1CorrelationValue[[i]],
                group2CorrelationValue = edgeData$group2CorrelationValue[[i]],
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

.computeDegreeMatchedScores <- function(rawRewiringScore, totalConnections,
                                        minimumBinSize = 5L) {
    logDegree <- log2(totalConnections + 1)
    degreeBins <- cut(
        x = logDegree,
        breaks = seq(from = 0, to = max(logDegree) + 1, by = 1),
        include.lowest = TRUE
    )
    zScores <- rep(NA_real_, length(rawRewiringScore))

    for (binName in levels(degreeBins)) {
        binIndex <- which(degreeBins == binName)
        binValues <- rawRewiringScore[binIndex]
        # Too few nodes cannot support a within-bin standardisation: with two
        # nodes scale() returns +/-0.707 whatever the data says. Leave those
        # missing rather than reporting a number that only looks like evidence.
        if (length(binIndex) < minimumBinSize || length(unique(binValues)) == 1L) {
            next
        }
        zScores[binIndex] <- as.numeric(scale(binValues))
    }

    list(
        degreeBin = degreeBins,
        degreeMatchedZScore = zScores
    )
}

.checkScoringEdges <- function(edgeData) {
    keyOne <- .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier)
    keyTwo <- .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier)
    if (any(keyOne == keyTwo)) {
        stop("Rewiring scores cannot be calculated from self-loops.", call. = FALSE)
    }

    edgeKey <- paste(pmin(keyOne, keyTwo), pmax(keyOne, keyTwo), sep = "\r")
    if (anyDuplicated(edgeKey)) {
        stop(
            "Rewiring scores require one row per edge; duplicate or mirrored edges were found.",
            call. = FALSE
        )
    }

    invisible(NULL)
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
#' @details These are descriptive CorNetto scores, not tests. None of them is
#'   a p-value, and `degreeMatchedZScore` is not a z-score against any null
#'   model. [permuteRewiringScores()] supplies a randomization p-value only
#'   under its documented fixed-universe and exchangeability conditions;
#'   otherwise it returns conditional rankings.
#'
#'   For a node with incident differential z-scores `z_1, ..., z_k`:
#'
#'   - `rawRewiringScore` is the L2 norm `sqrt(sum(z_i^2))`. It grows with
#'     degree: if the `z_i` were independent standard normals, `sum(z_i^2)`
#'     would follow a chi-square distribution on `k` degrees of freedom. Its
#'     squared score would then have expectation `k`, and the score itself
#'     would grow roughly as `sqrt(k)`. Edges sharing a node are not
#'     independent, so treat that only as intuition for the scaling.
#'   - `rootMeanSquareRewiringScore` divides by `sqrt(k)`, which removes most
#'     of that degree dependence. Under the same rough argument its squared
#'     value has expectation 1; the score itself is biased below 1 at small
#'     `k` and approaches 1 as degree grows. This is the most comparable score
#'     across nodes of different degree.
#'   - `degreeMatchedZScore` standardizes `rawRewiringScore` within bins of
#'     `log2(degree + 1)`. It is a within-bin ranking aid only. Roughly half
#'     the nodes in every bin are negative by construction, and the scale is
#'     not comparable across bins, so a bin holding few nodes cannot produce a
#'     meaningful spread. Bins with fewer than five nodes, or with no spread,
#'     return `NA`.
#'
#'   The network must contain one row per edge and no self-loops. Duplicate or
#'   mirrored edges would change both the degree and the score, so the function
#'   stops when either is detected. Score the network before calling
#'   [combineNetworks()] with `includeReverseEdges = TRUE`.
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
    .assertScalarLogical(storeResult, "storeResult")
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
    if (anyNA(edgeData$edgeWeight) || any(!is.finite(edgeData$edgeWeight))) {
        stop(
            "`differentialCorrelationNetwork$edgeWeight` must contain finite, non-missing values.",
            call. = FALSE
        )
    }
    .checkScoringEdges(edgeData)
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
    weightsByNode <- split(
        longTable$edgeWeight,
        factor(longTable$nodeKey, levels = nodeKeys)
    )
    firstRow <- match(nodeKeys, longTable$nodeKey)

    totalConnections <- lengths(weightsByNode)
    rawRewiringScore <- vapply(
        weightsByNode,
        function(w) sqrt(sum(w ^ 2)),
        numeric(1L)
    )

    rewiringTable <- data.frame(
        nodeKey = nodeKeys,
        nodeIdentifier = longTable$nodeIdentifier[firstRow],
        nodeName = longTable$nodeName[firstRow],
        assayName = longTable$assayName[firstRow],
        totalConnections = as.integer(totalConnections),
        rawRewiringScore = as.numeric(rawRewiringScore),
        rootMeanSquareRewiringScore = as.numeric(rawRewiringScore / sqrt(totalConnections)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )

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
