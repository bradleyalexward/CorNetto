#' Read an Assay Table from Disk
#'
#' Read a normalized abundance table from a delimited text file and
#' convert it to a `SummarizedExperiment`.
#'
#' @param filePath Path to a comma-separated or tab-separated text file.
#' @param assayName Optional assay label used only in messages.
#' @param featureIdentifierColumn Column containing unique feature
#'   identifiers. May be a column name or numeric index.
#' @param featureNameColumn Optional column containing human-readable
#'   feature labels. May be a column name or numeric index.
#' @param delimiter Optional delimiter. When `NULL`, the delimiter is
#'   inferred from the file extension.
#'
#' @return A `SummarizedExperiment`.
#' @export
readAssayData <- function(
    filePath,
    assayName = NULL,
    featureIdentifierColumn = 1L,
    featureNameColumn = NULL,
    delimiter = NULL
) {
    .assertScalarCharacter(filePath, "filePath")
    if (!file.exists(filePath)) {
        stop("`filePath` does not exist: ", filePath, call. = FALSE)
    }

    delimiter <- .matchDelimiter(filePath, delimiter)
    assayTable <- utils::read.table(
        file = filePath,
        header = TRUE,
        sep = delimiter,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (is.numeric(featureIdentifierColumn)) {
        featureIdentifierColumn <- names(assayTable)[featureIdentifierColumn]
    }

    if (!featureIdentifierColumn %in% names(assayTable)) {
        stop(
            "`featureIdentifierColumn` was not found in the input table.",
            call. = FALSE
        )
    }

    featureIdentifiers <- as.character(assayTable[[featureIdentifierColumn]])
    if (anyNA(featureIdentifiers) || any(!nzchar(featureIdentifiers))) {
        stop("Feature identifiers must be non-missing and non-empty.", call. = FALSE)
    }

    featureNames <- featureIdentifiers
    if (!is.null(featureNameColumn)) {
        if (is.numeric(featureNameColumn)) {
            featureNameColumn <- names(assayTable)[featureNameColumn]
        }
        if (!featureNameColumn %in% names(assayTable)) {
            stop("`featureNameColumn` was not found in the input table.", call. = FALSE)
        }
        featureNames <- as.character(assayTable[[featureNameColumn]])
        featureNames[is.na(featureNames) | !nzchar(featureNames)] <- featureIdentifiers[
            is.na(featureNames) | !nzchar(featureNames)
        ]
    }

    metadataColumns <- c(featureIdentifierColumn, featureNameColumn)
    metadataColumns <- metadataColumns[!is.na(metadataColumns) & nzchar(metadataColumns)]
    sampleColumns <- setdiff(names(assayTable), metadataColumns)
    if (!length(sampleColumns)) {
        stop("No sample columns were found in the assay table.", call. = FALSE)
    }

    assayMatrix <- as.matrix(assayTable[, sampleColumns, drop = FALSE])
    storage.mode(assayMatrix) <- "numeric"
    rownames(assayMatrix) <- featureIdentifiers

    rowData <- S4Vectors::DataFrame(
        featureIdentifier = featureIdentifiers,
        featureName = featureNames,
        check.names = FALSE
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = assayMatrix),
        rowData = rowData
    )
}

.coerceSingleAssay <- function(
    assayObject,
    assayName,
    featureIdentifierColumn = NULL,
    featureNameColumn = NULL
) {
    if (methods::is(assayObject, "SummarizedExperiment")) {
        rowDataFrame <- as.data.frame(SummarizedExperiment::rowData(assayObject))
        assayMatrix <- SummarizedExperiment::assay(assayObject)
        storage.mode(assayMatrix) <- "numeric"
        if (is.null(rownames(assayMatrix))) {
            stop(
                "Assay `", assayName, "` must have row names representing feature identifiers.",
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
                rowData = S4Vectors::DataFrame(rowDataFrame, check.names = FALSE)
            )
        )
    }

    assayTable <- as.data.frame(assayObject, check.names = FALSE, stringsAsFactors = FALSE)

    if (!is.null(featureIdentifierColumn)) {
        if (is.numeric(featureIdentifierColumn)) {
            featureIdentifierColumn <- names(assayTable)[featureIdentifierColumn]
        }
        if (!featureIdentifierColumn %in% names(assayTable)) {
            stop(
                "Feature identifier column `", featureIdentifierColumn,
                "` was not found for assay `", assayName, "`.",
                call. = FALSE
            )
        }
        rownames(assayTable) <- as.character(assayTable[[featureIdentifierColumn]])
    }

    if (is.null(rownames(assayTable))) {
        stop(
            "Assay `", assayName,
            "` must have row names or a `featureIdentifierColumn`.",
            call. = FALSE
        )
    }

    featureIdentifiers <- rownames(assayTable)
    featureNames <- featureIdentifiers

    if (!is.null(featureNameColumn)) {
        if (is.numeric(featureNameColumn)) {
            featureNameColumn <- names(assayTable)[featureNameColumn]
        }
        if (!featureNameColumn %in% names(assayTable)) {
            stop(
                "Feature name column `", featureNameColumn,
                "` was not found for assay `", assayName, "`.",
                call. = FALSE
            )
        }
        featureNames <- as.character(assayTable[[featureNameColumn]])
        featureNames[is.na(featureNames) | !nzchar(featureNames)] <- featureIdentifiers[
            is.na(featureNames) | !nzchar(featureNames)
        ]
    }

    dropColumns <- c(featureIdentifierColumn, featureNameColumn)
    dropColumns <- dropColumns[!is.na(dropColumns)]
    abundanceTable <- assayTable[, setdiff(names(assayTable), dropColumns), drop = FALSE]
    abundanceMatrix <- as.matrix(abundanceTable)
    storage.mode(abundanceMatrix) <- "numeric"
    rownames(abundanceMatrix) <- featureIdentifiers

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
#' @param featureIdentifierColumns Optional scalar or named list
#'   specifying feature identifier columns for data-frame assays.
#' @param featureNameColumns Optional scalar or named list specifying
#'   feature name columns for data-frame assays.
#'
#' @return A validated `MultiAssayExperiment`.
#' @export
createAnalysisData <- function(
    assayList,
    sampleData,
    featureIdentifierColumns = NULL,
    featureNameColumns = NULL
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
        featureIdentifierColumn <- featureIdentifierColumns
        featureNameColumn <- featureNameColumns

        if (is.list(featureIdentifierColumns) && !is.null(featureIdentifierColumns[[assayName]])) {
            featureIdentifierColumn <- featureIdentifierColumns[[assayName]]
        }
        if (is.list(featureNameColumns) && !is.null(featureNameColumns[[assayName]])) {
            featureNameColumn <- featureNameColumns[[assayName]]
        }

        experiments[[assayName]] <- .coerceSingleAssay(
            assayObject = assayList[[assayName]],
            assayName = assayName,
            featureIdentifierColumn = featureIdentifierColumn,
            featureNameColumn = featureNameColumn
        )
    }

    analysisData <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiments,
        colData = sampleData
    )

    analysisData <- .initializeCorNettoStore(analysisData)
    validateAnalysisData(analysisData)
}

#' Validate a CorNetto Analysis Object
#'
#' Validate that a `MultiAssayExperiment` is compatible with CorNetto.
#'
#' @param analysisData A `MultiAssayExperiment`.
#'
#' @return The validated `MultiAssayExperiment`.
#' @export
validateAnalysisData <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)

    sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData), stringsAsFactors = FALSE)
    if (is.null(rownames(sampleData)) || anyDuplicated(rownames(sampleData))) {
        stop("`colData` row names must be unique sample identifiers.", call. = FALSE)
    }

    experimentList <- MultiAssayExperiment::experiments(analysisData)
    if (!length(experimentList)) {
        stop("`analysisData` contains no assays.", call. = FALSE)
    }

    for (assayName in names(experimentList)) {
        assayObject <- experimentList[[assayName]]
        assayMatrix <- .extractAssayMatrix(analysisData, assayName)

        if (!is.numeric(assayMatrix)) {
            stop("Assay `", assayName, "` must contain numeric values.", call. = FALSE)
        }
        if (is.null(rownames(assayMatrix)) || anyDuplicated(rownames(assayMatrix))) {
            stop(
                "Assay `", assayName, "` must have unique feature identifiers as row names.",
                call. = FALSE
            )
        }
        if (is.null(colnames(assayMatrix)) || anyDuplicated(colnames(assayMatrix))) {
            stop(
                "Assay `", assayName, "` must have unique sample identifiers as column names.",
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
            rowDataFrame <- as.data.frame(SummarizedExperiment::rowData(assayObject))
            if (!"featureIdentifier" %in% names(rowDataFrame)) {
                rowDataFrame$featureIdentifier <- rownames(assayMatrix)
            }
            if (!"featureName" %in% names(rowDataFrame)) {
                rowDataFrame$featureName <- rowDataFrame$featureIdentifier
            }

            experimentList[[assayName]] <- SummarizedExperiment::SummarizedExperiment(
                assays = list(abundance = assayMatrix),
                rowData = S4Vectors::DataFrame(rowDataFrame, check.names = FALSE)
            )
        }
    }

    MultiAssayExperiment::experiments(analysisData) <- experimentList
    .initializeCorNettoStore(analysisData)
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
#'   minimum variance required for retention.
#'
#' @return A filtered `MultiAssayExperiment`.
#' @export
filterFeatures <- function(
    analysisData,
    assayNames = NULL,
    featureSubset = NULL,
    topVariableFeatures = NULL,
    minimumVariance = NULL
) {
    analysisData <- validateAnalysisData(analysisData)
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
            keepFeatures <- intersect(keepFeatures, as.character(assayFeatureSubset))
        }

        assayMinimumVariance <- minimumVariance
        if (is.list(minimumVariance) && !is.null(minimumVariance[[assayName]])) {
            assayMinimumVariance <- minimumVariance[[assayName]]
        }
        if (!is.null(assayMinimumVariance)) {
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
    validateAnalysisData(analysisData)
}

#' Map Feature Identifiers to Display Names
#'
#' Retrieve assay-specific feature identifiers and their labels.
#'
#' @param analysisData A `MultiAssayExperiment`.
#' @param assayName Optional assay name. When `NULL`, mappings are
#'   returned for all assays.
#' @param featureIdentifiers Optional feature identifiers to subset.
#' @param featureNameColumn Row-data column containing display names.
#'
#' @return A `DataFrame` with assay-specific feature mappings.
#' @export
mapFeatureIdentifiers <- function(
    analysisData,
    assayName = NULL,
    featureIdentifiers = NULL,
    featureNameColumn = "featureName"
) {
    analysisData <- validateAnalysisData(analysisData)
    assayNames <- assayName
    if (is.null(assayNames)) {
        assayNames <- names(MultiAssayExperiment::experiments(analysisData))
    }

    mappingList <- lapply(
        assayNames,
        function(singleAssayName) {
            .extractFeatureAnnotations(
                analysisData = analysisData,
                assayName = singleAssayName,
                featureNameColumn = featureNameColumn
            )
        }
    )

    mappingTable <- do.call(rbind, mappingList)
    if (!is.null(featureIdentifiers)) {
        mappingTable <- mappingTable[
            mappingTable$featureIdentifier %in% as.character(featureIdentifiers),
            ,
            drop = FALSE
        ]
    }

    S4Vectors::DataFrame(mappingTable, check.names = FALSE)
}
