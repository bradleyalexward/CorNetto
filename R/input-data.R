.coerceSingleAssay <- function(
    assayObject,
    assayName
) {
    if (methods::is(assayObject, "SummarizedExperiment")) {
        rowDataFrame <- .asPlainDataFrame(
            SummarizedExperiment::rowData(assayObject)
        )
        assayMatrix <- SummarizedExperiment::assay(assayObject)
        assayMatrix <- .coerceNumericMatrixStrict(
            assayMatrix,
            context = paste0("assay `", assayName, "`")
        )
        if (is.null(rownames(assayMatrix))) {
            stop(
                "Assay `", assayName,
                "` must have row names representing feature identifiers.",
                call. = FALSE
            )
        }

        if (!"featureIdentifier" %in% names(rowDataFrame)) {
            rowDataFrame$featureIdentifier <- rownames(assayMatrix)
        }
        if (!"featureName" %in% names(rowDataFrame)) {
            rowDataFrame$featureName <- rowDataFrame$featureIdentifier
        }

        return(
            SummarizedExperiment::SummarizedExperiment(
                assays = list(abundance = assayMatrix),
                rowData = S4Vectors::DataFrame(
                    rowDataFrame,
                    check.names = FALSE
                )
            )
        )
    }

    assayTable <- .asPlainDataFrame(assayObject)

    if (.hasAutomaticRowNames(assayTable)) {
        stop(
            "Assay `", assayName,
            "` must have row names representing feature identifiers.",
            call. = FALSE
        )
    }

    featureIdentifiers <- rownames(assayTable)
    if (anyNA(featureIdentifiers) || any(!nzchar(featureIdentifiers)) ||
        anyDuplicated(featureIdentifiers)) {
        stop(
            "Assay `", assayName,
            "` must have unique, non-missing feature identifiers.",
            call. = FALSE
        )
    }
    featureNames <- featureIdentifiers

    abundanceTable <- assayTable
    abundanceMatrix <- as.matrix(abundanceTable)
    rownames(abundanceMatrix) <- featureIdentifiers
    abundanceMatrix <- .coerceNumericMatrixStrict(
        abundanceMatrix,
        context = paste0("assay `", assayName, "`")
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = abundanceMatrix),
        rowData = S4Vectors::DataFrame(
            featureIdentifier = featureIdentifiers,
            featureName = featureNames,
            check.names = FALSE
        )
    )
}

#' Create a CorNetto Analysis Object
#'
#' Create a `MultiAssayExperiment` from normalized assay objects and
#' sample metadata.
#'
#' @param assayList A named list of assays. Each element may be a matrix,
#'   data frame, or `SummarizedExperiment`.
#' @param sampleData Sample metadata with row names corresponding to
#'   sample identifiers. A `sampleId` column may be used when row names
#'   are absent.
#'
#' @return A validated `MultiAssayExperiment`.
#' @examples
#' protein <- matrix(
#'     seq_len(12),
#'     nrow = 3,
#'     dimnames = list(paste0("P", 1:3), paste0("S", 1:4))
#' )
#' sampleData <- data.frame(
#'     sampleId = paste0("S", 1:4),
#'     group = rep(c("A", "B"), each = 2)
#' )
#' analysisData <- createAnalysisData(list(protein = protein), sampleData)
#' analysisData
#' @export
createAnalysisData <- function(
    assayList,
    sampleData
) {
    if (!is.list(assayList) || !length(assayList)) {
        stop("`assayList` must be a named list of assays.", call. = FALSE)
    }
    if (is.null(names(assayList)) || any(!nzchar(names(assayList)))) {
        stop("`assayList` must have non-empty assay names.", call. = FALSE)
    }

    sampleData <- .coerceSampleData(sampleData)
    experiments <- vector("list", length(assayList))
    names(experiments) <- names(assayList)

    for (assayName in names(assayList)) {
        experiments[[assayName]] <- .coerceSingleAssay(
            assayObject = assayList[[assayName]],
            assayName = assayName
        )
    }

    analysisData <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiments,
        colData = sampleData
    )

    analysisData <- .initializeCorNettoStore(analysisData)
    validateAnalysisData(analysisData, quiet = TRUE)
}

#' Validate a CorNetto Analysis Object
#'
#' Validate that a `MultiAssayExperiment` is compatible with CorNetto.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param quiet Whether to suppress the success message.
#'
#' @return The validated `MultiAssayExperiment`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' validated <- validateAnalysisData(analysisData)
#' validated
#' @export
validateAnalysisData <- function(analysisData, quiet = FALSE) {
    .assertMultiAssayExperiment(analysisData)
    .assertScalarLogical(quiet, "quiet")

    sampleData <- .asPlainDataFrame(
        MultiAssayExperiment::colData(analysisData)
    )
    if (is.null(rownames(sampleData)) || anyNA(rownames(sampleData)) ||
        any(!nzchar(rownames(sampleData))) ||
        anyDuplicated(rownames(sampleData))) {
        stop(
            "`colData` row names must be unique, non-missing sample ",
            "identifiers.",
            call. = FALSE
        )
    }

    experimentList <- MultiAssayExperiment::experiments(analysisData)
    if (!length(experimentList)) {
        stop("`analysisData` contains no assays.", call. = FALSE)
    }

    for (assayName in names(experimentList)) {
        assayObject <- experimentList[[assayName]]
        assayMatrix <- .extractAssayMatrix(analysisData, assayName)

        if (!is.numeric(assayMatrix)) {
            stop(
                "Assay `", assayName, "` must contain numeric values.",
                call. = FALSE
            )
        }
        if (is.null(rownames(assayMatrix)) ||
            anyNA(rownames(assayMatrix)) ||
            any(!nzchar(rownames(assayMatrix))) ||
            anyDuplicated(rownames(assayMatrix))) {
            stop(
                "Assay `", assayName,
                "` must have unique, non-missing feature identifiers ",
                "as row names.",
                call. = FALSE
            )
        }
        if (is.null(colnames(assayMatrix)) ||
            anyNA(colnames(assayMatrix)) ||
            any(!nzchar(colnames(assayMatrix))) ||
            anyDuplicated(colnames(assayMatrix))) {
            stop(
                "Assay `", assayName,
                "` must have unique, non-missing sample identifiers ",
                "as column names.",
                call. = FALSE
            )
        }
        if (!all(colnames(assayMatrix) %in% rownames(sampleData))) {
            stop(
                "All assay samples in `", assayName,
                "` must be present in `colData` row names.",
                call. = FALSE
            )
        }

        if (methods::is(assayObject, "SummarizedExperiment")) {
            rowDataFrame <- .asPlainDataFrame(
                SummarizedExperiment::rowData(assayObject)
            )
            if (!"featureIdentifier" %in% names(rowDataFrame)) {
                rowDataFrame$featureIdentifier <- rownames(assayMatrix)
            }
            if (!"featureName" %in% names(rowDataFrame)) {
                rowDataFrame$featureName <- rowDataFrame$featureIdentifier
            }

            experimentList[[assayName]] <-
                SummarizedExperiment::SummarizedExperiment(
                    assays = list(abundance = assayMatrix),
                    rowData = S4Vectors::DataFrame(
                        rowDataFrame,
                        check.names = FALSE
                    )
                )
        }
    }

    MultiAssayExperiment::experiments(analysisData) <- experimentList
    analysisData <- .initializeCorNettoStore(analysisData)
    if (!isTRUE(quiet)) {
        message("Analysis data object is structured correctly for CorNetto.")
    }
    analysisData
}

#' Summarize a CorNetto Analysis Object
#'
#' Summarize assay dimensions, missingness, low-variance features, sample
#' overlap, and optional group sizes before network construction.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param groupColumn Optional sample metadata column used to summarize
#'   group sizes.
#' @param missingnessWarningThreshold Fraction of missing assay values
#'   above which a warning is reported for an assay.
#' @param lowSampleWarningThreshold Group or assay sample count below
#'   which a warning is reported.
#' @param highFeatureWarningThreshold Feature count above which dense
#'   all-pairs correlation is flagged as computationally heavy.
#' @param quiet Whether to suppress emitted warnings and only return them
#'   in the result object.
#'
#' @return A named list with `assays`, `samples`, and `warnings`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' summary <- summarizeAnalysisData(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     quiet = TRUE
#' )
#' summary$assays
#' @export
summarizeAnalysisData <- function(
    analysisData,
    groupColumn = NULL,
    missingnessWarningThreshold = 0.2,
    lowSampleWarningThreshold = 10L,
    highFeatureWarningThreshold = 3000L,
    quiet = FALSE
) {
    .assertMultiAssayExperiment(analysisData)
    .assertScalarLogical(quiet, "quiet")

    sampleData <- .asPlainDataFrame(MultiAssayExperiment::colData(analysisData))
    experimentList <- MultiAssayExperiment::experiments(analysisData)
    warningMessages <- character()

    assayResults <- lapply(
        names(experimentList),
        function(assayName) {
            assayWarnings <- character()
            assayMatrix <- .extractAssayMatrix(analysisData, assayName)
            featureCount <- nrow(assayMatrix)
            sampleCount <- ncol(assayMatrix)
            missingCount <- sum(is.na(assayMatrix))
            totalValues <- length(assayMatrix)
            missingFraction <- if (totalValues) missingCount / totalValues else NA_real_
            featureVariance <- apply(assayMatrix, 1L, stats::var, na.rm = TRUE)
            allMissingFeatureCount <- sum(rowSums(!is.na(assayMatrix)) == 0L)
            zeroVarianceFeatureCount <- sum(!is.na(featureVariance) & featureVariance == 0)
            highMissingFeatureCount <- sum(rowMeans(is.na(assayMatrix)) > missingnessWarningThreshold)
            samplesMissingFromColData <- setdiff(colnames(assayMatrix), rownames(sampleData))

            if (sampleCount < lowSampleWarningThreshold) {
                assayWarnings <- c(
                    assayWarnings,
                    paste0(
                        "Assay `", assayName, "` has ", sampleCount,
                        " samples; correlation estimates are unstable below ",
                        lowSampleWarningThreshold, " samples."
                    )
                )
            }
            if (featureCount > highFeatureWarningThreshold) {
                assayWarnings <- c(
                    assayWarnings,
                    paste0(
                        "Assay `", assayName, "` has ", featureCount,
                        " features; dense all-pairs correlation can be slow and memory intensive."
                    )
                )
            }
            if (!is.na(missingFraction) && missingFraction > missingnessWarningThreshold) {
                assayWarnings <- c(
                    assayWarnings,
                    paste0(
                        "Assay `", assayName, "` has ",
                        round(100 * missingFraction, 1),
                        "% missing values."
                    )
                )
            }
            if (zeroVarianceFeatureCount > 0L) {
                assayWarnings <- c(
                    assayWarnings,
                    paste0(
                        "Assay `", assayName, "` has ", zeroVarianceFeatureCount,
                        " zero-variance features."
                    )
                )
            }
            if (length(samplesMissingFromColData)) {
                assayWarnings <- c(
                    assayWarnings,
                    paste0(
                        "Assay `", assayName, "` has samples missing from `colData`: ",
                        paste(utils::head(samplesMissingFromColData, 5L), collapse = ", ")
                    )
                )
            }

            list(
                row = data.frame(
                    assayName = assayName,
                    featureCount = featureCount,
                    sampleCount = sampleCount,
                    missingValueCount = missingCount,
                    missingValueFraction = missingFraction,
                    zeroVarianceFeatureCount = zeroVarianceFeatureCount,
                    allMissingFeatureCount = allMissingFeatureCount,
                    highMissingFeatureCount = highMissingFeatureCount,
                    samplesMissingFromColData = length(samplesMissingFromColData),
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                ),
                warnings = assayWarnings
            )
        }
    )
    assayRows <- lapply(assayResults, `[[`, "row")
    warningMessages <- unlist(lapply(assayResults, `[[`, "warnings"),
        use.names = FALSE)
    assaySummary <- do.call(rbind, assayRows)

    assaySamples <- unique(unlist(
        lapply(
            names(experimentList),
            function(assayName) colnames(.extractAssayMatrix(analysisData, assayName))
        ),
        use.names = FALSE
    ))
    sampleSummary <- data.frame(
        sampleCount = nrow(sampleData),
        samplesInAnyAssay = sum(rownames(sampleData) %in% assaySamples),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!is.null(groupColumn)) {
        .assertScalarCharacter(groupColumn, "groupColumn")
        if (!groupColumn %in% names(sampleData)) {
            stop("`groupColumn` was not found in `colData`: ", groupColumn, call. = FALSE)
        }
        groupCounts <- as.data.frame(table(sampleData[[groupColumn]]), stringsAsFactors = FALSE)
        names(groupCounts) <- c("groupLevel", "sampleCount")
        lowGroups <- groupCounts$groupLevel[groupCounts$sampleCount < lowSampleWarningThreshold]
        if (length(lowGroups)) {
            warningMessages <- c(
                warningMessages,
                paste0(
                    "The following groups in `", groupColumn, "` have fewer than ",
                    lowSampleWarningThreshold, " samples: ",
                    paste(lowGroups, collapse = ", ")
                )
            )
        }
        sampleSummary <- list(
            overall = S4Vectors::DataFrame(sampleSummary, check.names = FALSE),
            groups = S4Vectors::DataFrame(groupCounts, check.names = FALSE)
        )
    } else {
        sampleSummary <- list(
            overall = S4Vectors::DataFrame(sampleSummary, check.names = FALSE),
            groups = NULL
        )
    }

    warningMessages <- unique(warningMessages)
    if (!isTRUE(quiet) && length(warningMessages)) {
        warning(paste(warningMessages, collapse = "\n"), call. = FALSE)
    }

    list(
        assays = S4Vectors::DataFrame(assaySummary, check.names = FALSE),
        samples = sampleSummary,
        warnings = warningMessages
    )
}

#' Filter Samples Before Network Construction
#'
#' Filter a CorNetto analysis object to selected samples while allowing
#' assays to retain different overlapping subsets of those samples.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param sampleIds Optional sample identifiers to retain.
#' @param groupColumn Optional sample metadata column used to retain
#'   samples by group membership.
#' @param groupLevels Optional group labels to retain from `groupColumn`.
#'   When `NULL` and `groupColumn` is supplied, samples with non-missing
#'   and non-empty values in `groupColumn` are retained.
#' @param dropUnusedLevels Whether to drop unused factor levels in
#'   sample metadata after filtering.
#'
#' @return A filtered `MultiAssayExperiment`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' filtered <- filterSamples(
#'     analysisData,
#'     groupColumn = "clinicalGroup",
#'     groupLevels = "PASC"
#' )
#' summarizeAnalysisData(
#'     filtered,
#'     groupColumn = "clinicalGroup",
#'     quiet = TRUE
#' )$samples
#' @export
filterSamples <- function(
    analysisData,
    sampleIds = NULL,
    groupColumn = NULL,
    groupLevels = NULL,
    dropUnusedLevels = TRUE
) {
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    .assertScalarLogical(dropUnusedLevels, "dropUnusedLevels")

    if (is.null(sampleIds) && is.null(groupColumn)) {
        stop(
            "Supply `sampleIds`, `groupColumn`, or both to filter samples.",
            call. = FALSE
        )
    }

    sampleData <- .asPlainDataFrame(MultiAssayExperiment::colData(analysisData))
    selectedSampleIds <- rownames(sampleData)

    if (!is.null(sampleIds)) {
        sampleIds <- as.character(sampleIds)
        if (anyNA(sampleIds) || any(!nzchar(sampleIds))) {
            stop(
                "`sampleIds` must be non-missing and non-empty.",
                call. = FALSE
            )
        }
        missingSampleIds <- setdiff(sampleIds, rownames(sampleData))
        if (length(missingSampleIds)) {
            stop(
                "`sampleIds` were not found in `colData`: ",
                paste(utils::head(missingSampleIds, 5L), collapse = ", "),
                call. = FALSE
            )
        }
        selectedSampleIds <- intersect(selectedSampleIds, sampleIds)
    }

    if (!is.null(groupColumn)) {
        .assertScalarCharacter(groupColumn, "groupColumn")
        if (!groupColumn %in% names(sampleData)) {
            stop(
                "`groupColumn` was not found in `colData`: ",
                groupColumn,
                call. = FALSE
            )
        }

        groupValues <- as.character(sampleData[[groupColumn]])
        if (is.null(groupLevels)) {
            groupIndex <- !is.na(groupValues) & nzchar(groupValues)
        } else {
            groupLevels <- as.character(groupLevels)
            badGroupLevels <- !length(groupLevels) ||
                anyNA(groupLevels) ||
                any(!nzchar(groupLevels))
            if (badGroupLevels) {
                stop(
                    "`groupLevels` must be non-missing and non-empty.",
                    call. = FALSE
                )
            }
            groupIndex <- groupValues %in% groupLevels
        }
        selectedSampleIds <- intersect(
            selectedSampleIds,
            rownames(sampleData)[groupIndex]
        )
    }

    if (!length(selectedSampleIds)) {
        stop("Sample filtering removed all samples.", call. = FALSE)
    }

    analysisData <- analysisData[, selectedSampleIds, drop = FALSE]

    sampleData <- .asPlainDataFrame(MultiAssayExperiment::colData(analysisData))
    sampleData <- sampleData[selectedSampleIds, , drop = FALSE]
    if (isTRUE(dropUnusedLevels)) {
        sampleData[] <- lapply(
            sampleData,
            function(sampleColumn) {
                if (is.factor(sampleColumn)) {
                    return(droplevels(sampleColumn))
                }
                sampleColumn
            }
        )
    }
    sampleData <- S4Vectors::DataFrame(sampleData, check.names = FALSE)
    rownames(sampleData) <- selectedSampleIds

    MultiAssayExperiment::colData(analysisData) <- sampleData
    validateAnalysisData(analysisData, quiet = TRUE)
}

#' Filter Features Before Network Construction
#'
#' Filter one or more assays to reduce network complexity.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param assayNames Optional assay names to filter. When `NULL`, all
#'   assays are filtered.
#' @param featureSubset Optional vector or named list of feature
#'   identifiers to retain.
#' @param topVariableFeatures Optional scalar or named list giving the
#'   number of most variable features to retain per assay.
#' @param minimumVariance Optional scalar or named list defining the
#'   minimum variance required for retention. Features are retained when
#'   their variance is greater than or equal to this value.
#' @param removeZeroVariance Logical scalar or named list. When `TRUE`,
#'   remove features with zero or undefined variance before applying
#'   `minimumVariance` or `topVariableFeatures`.
#'
#' @return A filtered `MultiAssayExperiment`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' filtered <- filterFeatures(
#'     analysisData,
#'     assayNames = "protein",
#'     topVariableFeatures = 3
#' )
#' filtered
#' @export
filterFeatures <- function(
    analysisData,
    assayNames = NULL,
    featureSubset = NULL,
    topVariableFeatures = NULL,
    minimumVariance = NULL,
    removeZeroVariance = FALSE
) {
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
    experimentList <- MultiAssayExperiment::experiments(analysisData)

    if (is.null(assayNames)) {
        assayNames <- names(experimentList)
    }

    for (assayName in assayNames) {
        assayObject <- experimentList[[assayName]]
        assayMatrix <- .extractAssayMatrix(analysisData, assayName)
        keepFeatures <- rownames(assayMatrix)

        assayFeatureSubset <- featureSubset
        if (is.list(featureSubset) && !is.null(featureSubset[[assayName]])) {
            assayFeatureSubset <- featureSubset[[assayName]]
        }
        if (!is.null(assayFeatureSubset)) {
            keepFeatures <- intersect(
                keepFeatures,
                as.character(assayFeatureSubset)
            )
        }

        assayRemoveZeroVariance <- removeZeroVariance
        if (is.list(removeZeroVariance)) {
            assayRemoveZeroVariance <- FALSE
            if (!is.null(removeZeroVariance[[assayName]])) {
                assayRemoveZeroVariance <- removeZeroVariance[[assayName]]
            }
        }
        .assertScalarLogical(assayRemoveZeroVariance, "removeZeroVariance")
        if (isTRUE(assayRemoveZeroVariance)) {
            featureVariance <- apply(assayMatrix, 1L, stats::var, na.rm = TRUE)
            keepFeatures <- intersect(
                keepFeatures,
                names(featureVariance)[
                    !is.na(featureVariance) & featureVariance > 0
                ]
            )
        }

        assayMinimumVariance <- minimumVariance
        if (is.list(minimumVariance) &&
            !is.null(minimumVariance[[assayName]])) {
            assayMinimumVariance <- minimumVariance[[assayName]]
        }
        if (!is.null(assayMinimumVariance)) {
            if (!is.numeric(assayMinimumVariance) ||
                length(assayMinimumVariance) != 1L ||
                is.na(assayMinimumVariance) || assayMinimumVariance < 0) {
                stop(
                    "`minimumVariance` must be a non-negative numeric scalar.",
                    call. = FALSE
                )
            }
            featureVariance <- apply(assayMatrix, 1L, stats::var, na.rm = TRUE)
            keepFeatures <- intersect(
                keepFeatures,
                names(featureVariance)[featureVariance >= assayMinimumVariance]
            )
        }

        assayTopVariableFeatures <- topVariableFeatures
        if (is.list(topVariableFeatures) && !is.null(topVariableFeatures[[assayName]])) {
            assayTopVariableFeatures <- topVariableFeatures[[assayName]]
        }
        if (!is.null(assayTopVariableFeatures)) {
            if (!is.numeric(assayTopVariableFeatures) || length(assayTopVariableFeatures) != 1L ||
                is.na(assayTopVariableFeatures) || assayTopVariableFeatures < 1L ||
                assayTopVariableFeatures != as.integer(assayTopVariableFeatures)) {
                stop("`topVariableFeatures` must be a positive integer.", call. = FALSE)
            }
            featureVariance <- apply(
                assayMatrix[keepFeatures, , drop = FALSE],
                1L,
                stats::var,
                na.rm = TRUE
            )
            featureVariance <- sort(featureVariance, decreasing = TRUE)
            keepFeatures <- names(utils::head(featureVariance, assayTopVariableFeatures))
        }

        keepFeatures <- intersect(rownames(assayMatrix), keepFeatures)
        if (!length(keepFeatures)) {
            stop("Filtering removed all features from assay `", assayName, "`.", call. = FALSE)
        }

        experimentList[[assayName]] <- assayObject[keepFeatures, , drop = FALSE]
    }

    MultiAssayExperiment::experiments(analysisData) <- experimentList
    validateAnalysisData(analysisData, quiet = TRUE)
}
