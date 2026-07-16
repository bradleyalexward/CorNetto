.computePearsonCorrelations <- function(assayMatrix) {
    if (nrow(assayMatrix) < 2L) {
        return(list(
            correlationMatrix = matrix(numeric(), nrow = 0L, ncol = 0L),
            pValueMatrix = matrix(numeric(), nrow = 0L, ncol = 0L),
            sampleCountMatrix = matrix(integer(), nrow = 0L, ncol = 0L)
        ))
    }

    correlationMatrix <- stats::cor(
        x = t(assayMatrix),
        use = "pairwise.complete.obs",
        method = "pearson"
    )
    sampleCountMatrix <- tcrossprod(!is.na(assayMatrix))
    statisticMatrix <- matrix(
        NA_real_,
        nrow = nrow(correlationMatrix),
        ncol = ncol(correlationMatrix),
        dimnames = dimnames(correlationMatrix)
    )
    testable <- !is.na(correlationMatrix) & sampleCountMatrix >= 3L
    perfect <- testable & abs(correlationMatrix) >= 1
    regular <- testable & !perfect
    statisticMatrix[perfect] <- sign(correlationMatrix[perfect]) * Inf
    statisticMatrix[regular] <- correlationMatrix[regular] * sqrt(
        (sampleCountMatrix[regular] - 2) /
            (1 - correlationMatrix[regular] ^ 2)
    )
    pValueMatrix <- matrix(
        NA_real_,
        nrow = nrow(correlationMatrix),
        ncol = ncol(correlationMatrix),
        dimnames = dimnames(correlationMatrix)
    )
    pValueMatrix[testable] <- 2 * stats::pt(
        q = abs(statisticMatrix[testable]),
        df = sampleCountMatrix[testable] - 2,
        lower.tail = FALSE
    )

    # cor() still returns +/-1 through two points, but there is no test to go
    # with it; drop both so an untestable pair cannot reach the edge table.
    correlationMatrix[sampleCountMatrix < 3] <- NA_real_
    pValueMatrix[sampleCountMatrix < 3] <- NA_real_
    diag(pValueMatrix) <- 0
    diag(correlationMatrix) <- 1

    list(
        correlationMatrix = correlationMatrix,
        pValueMatrix = pValueMatrix,
        sampleCountMatrix = sampleCountMatrix
    )
}

.computeCorrelationByPairwiseTests <- function(assayMatrix, correlationMethod) {
    featureCount <- nrow(assayMatrix)
    featureNames <- rownames(assayMatrix)
    if (featureCount < 2L) {
        return(list(
            correlationMatrix = matrix(numeric(), nrow = 0L, ncol = 0L),
            pValueMatrix = matrix(numeric(), nrow = 0L, ncol = 0L),
            sampleCountMatrix = matrix(integer(), nrow = 0L, ncol = 0L)
        ))
    }
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
    sampleCountMatrix <- matrix(
        NA_integer_,
        nrow = featureCount,
        ncol = featureCount,
        dimnames = list(featureNames, featureNames)
    )

    diag(correlationMatrix) <- 1
    diag(pValueMatrix) <- 0
    diag(sampleCountMatrix) <- rowSums(!is.na(assayMatrix))

    for (i in seq_len(featureCount - 1L)) {
        for (j in seq.int(i + 1L, featureCount)) {
            testResult <- .runCorrelationTest(
                x = as.numeric(assayMatrix[i, ]),
                y = as.numeric(assayMatrix[j, ]),
                correlationMethod = correlationMethod
            )
            if (is.null(testResult)) {
                next
            }

            correlationMatrix[i, j] <- testResult$correlationValue
            correlationMatrix[j, i] <- testResult$correlationValue
            pValueMatrix[i, j] <- testResult$pValue
            pValueMatrix[j, i] <- testResult$pValue
            sampleCountMatrix[i, j] <- testResult$sampleCount
            sampleCountMatrix[j, i] <- testResult$sampleCount
        }
    }

    list(
        correlationMatrix = correlationMatrix,
        pValueMatrix = pValueMatrix,
        sampleCountMatrix = sampleCountMatrix
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
    sampleCountMatrix = NULL,
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
    sampleCounts <- if (is.null(sampleCountMatrix)) {
        rep(NA_integer_, length(correlationValues))
    } else {
        as.integer(sampleCountMatrix[upperIndex])
    }
    adjustedPValues <- .adjustPValuesWithMissing(pValues, pAdjustMethod)

    edgeTable <- data.frame(
        fromFeatureIdentifier = fromIdentifiers,
        toFeatureIdentifier = toIdentifiers,
        correlationValue = correlationValues,
        pValue = pValues,
        adjustedPValue = adjustedPValues,
        sampleCount = sampleCounts,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    keepIndex <- !is.na(edgeTable$correlationValue)
    if (!is.null(minimumAbsoluteCorrelation)) {
        keepIndex <- keepIndex & abs(edgeTable$correlationValue) >= minimumAbsoluteCorrelation
    }
    if (!is.null(adjustedPValueThreshold)) {
        # An NA here means the pair was never testable, not that it failed the
        # threshold. Without the is.na() guard keepIndex carries NA and
        # subsetting a data frame with NA appends a row of all-missing values.
        keepIndex <- keepIndex &
            !is.na(edgeTable$adjustedPValue) &
            edgeTable$adjustedPValue <= adjustedPValueThreshold
    }
    edgeTable <- edgeTable[which(keepIndex), , drop = FALSE]

    if (!nrow(edgeTable)) {
        return(.emptyStandardEdgeTable())
    }

    featureAnnotations <- .asPlainDataFrame(featureAnnotations)
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
    edgeTable$group1SampleCount <- NA_real_
    edgeTable$group2SampleCount <- NA_real_
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
    if (length(sampleIds) < 10L) {
        warning(
            "Only ", length(sampleIds), " samples are available for assay `",
            assayName, "`. Correlation estimates are unstable below 10 samples.",
            call. = FALSE
        )
    }

    subsetMatrix <- assayMatrix[featureSubset, sampleIds, drop = FALSE]
    if (nrow(subsetMatrix) > 3000L) {
        warning(
            "Assay `", assayName, "` has ", nrow(subsetMatrix),
            " features selected for dense all-pairs correlation. Dense mode is intended ",
            "for small-to-medium assays and can be slow and memory intensive above 3000 features.",
            call. = FALSE
        )
    }

    subsetMatrix
}

#' Create a Group-Specific Correlation Network
#'
#' Create a within-omic correlation network for one assay and one sample
#' group. All pairs in the assay are tested, which is intended for
#' small-to-medium assays. For large assays, restrict the feature set with
#' [filterFeatures()] or `featureSubset`, or work from prior-supported pairs
#' with `testDifferentialCorrelation(candidateEdgeTable = ...)`.
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
#' @param pAdjustMethod Multiple-testing correction method. Use one of
#'   `stats::p.adjust.methods` or `"qvalue"`; the latter requires the suggested
#'   `qvalue` package.
#' @param featureNameColumn Row-data column containing display names.
#' @param resultName Optional name used when storing the result.
#' @param storeResult Whether to store the result in
#'   `metadata(analysisData)$cornetto`.
#'
#' @return A standardized edge `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @details Pearson p-values come from the t statistic
#'   \eqn{t = r * sqrt((n - 2) / (1 - r^2))} on \eqn{n - 2} degrees of
#'   freedom, matching
#'   [stats::cor.test()]. Spearman and Kendall p-values are taken from
#'   [stats::cor.test()] with `exact = FALSE`, using its asymptotic
#'   approximations. Tied ranks are permitted; Kendall's estimate is tau-b.
#'
#'   Correlations use pairwise complete observations, so `n` varies by pair and
#'   is reported per edge in `sampleCount`. Pairs with fewer than three shared
#'   observations, and pairs where either feature is constant, have no test and
#'   are dropped rather than reported as non-significant.
#'
#'   The Pearson t test assumes independent observations and approximate
#'   bivariate normality for each pair. Pairwise deletion is defensible only
#'   when the missingness process is ignorable for the correlation of interest;
#'   abundance-dependent missingness in proteomics or metabolomics can violate
#'   that condition. Rank correlations reduce sensitivity to outliers but do
#'   not remove the independence or missingness assumptions.
#' @references
#' Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate.
#' \emph{Journal of the Royal Statistical Society B}, 57, 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Storey JD, Tibshirani R (2003). Statistical significance for genomewide
#' studies. \emph{PNAS}, 100, 9440-9445. \doi{10.1073/pnas.1530509100}
#' @examples
#' analysisData <- exampleAnalysisData()
#' correlationNetwork <- createCorrelationNetwork(
#'     analysisData,
#'     assayName = "protein",
#'     minimumAbsoluteCorrelation = 0,
#'     adjustedPValueThreshold = 1,
#'     storeResult = FALSE
#' )
#' correlationNetwork
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
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    correlationMethod <- match.arg(correlationMethod)
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
            slotName = "correlationResults",
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
        sampleCountMatrix = correlationResult$sampleCountMatrix,
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
        slotName = "correlationResults",
        resultObject = edgeTable,
        resultName = resultName
    )
}
