.corNettoStoreTemplate <- function() {
    list(
        knowledgeNetworks = list(),
        correlationResults = list(),
        differentialCorrelationResults = list(),
        differentialCorrelationNetworks = list(),
        rewiringResults = list(),
        validationResults = list(),
        integratedNetworks = list()
    )
}

.standardEdgeColumns <- function() {
    c(
        "fromFeatureIdentifier",
        "toFeatureIdentifier",
        "fromFeatureName",
        "toFeatureName",
        "fromAssayName",
        "toAssayName",
        "edgeType",
        "edgeDirection",
        "sourceType",
        "correlationScope",
        "correlationMethod",
        "knowledgeSource",
        "groupName",
        "comparisonName",
        "correlationValue",
        "group1CorrelationValue",
        "group2CorrelationValue",
        "pValue",
        "adjustedPValue",
        "group1PValue",
        "group2PValue",
        "group1AdjustedPValue",
        "group2AdjustedPValue",
        "zScoreDifference",
        "sampleCount",
        "group1SampleCount",
        "group2SampleCount",
        "edgeWeight",
        "evidenceScore",
        "isDirected"
    )
}

.numericEdgeColumns <- function() {
    c(
        "correlationValue",
        "group1CorrelationValue",
        "group2CorrelationValue",
        "pValue",
        "adjustedPValue",
        "group1PValue",
        "group2PValue",
        "group1AdjustedPValue",
        "group2AdjustedPValue",
        "zScoreDifference",
        "sampleCount",
        "group1SampleCount",
        "group2SampleCount",
        "edgeWeight",
        "evidenceScore"
    )
}

.logicalEdgeColumns <- function() {
    "isDirected"
}

## Adjust only the p-values that exist, and put NA back where there was no
## test to begin with. p.adjust() would otherwise count untestable pairs in
## the denominator.
.adjustPValuesWithMissing <- function(pValues, pAdjustMethod) {
    adjustedValues <- rep(NA_real_, length(pValues))
    keepIndex <- !is.na(pValues)
    if (!any(keepIndex)) {
        return(adjustedValues)
    }

    observedPValues <- as.numeric(pValues[keepIndex])
    if (identical(pAdjustMethod, "qvalue")) {
        if (!requireNamespace("qvalue", quietly = TRUE)) {
            stop(
                "The `qvalue` package is required when `pAdjustMethod = \"qvalue\"`.",
                call. = FALSE
            )
        }
        if (length(observedPValues) == 1L) {
            adjustedValues[keepIndex] <- observedPValues
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

        adjustedValues[keepIndex] <- qvalueResult$qvalues
        return(adjustedValues)
    }

    adjustedValues[keepIndex] <- stats::p.adjust(observedPValues, method = pAdjustMethod)
    adjustedValues
}

## One string per edge row, for exact-match comparisons between tables.
.edgeRowKey <- function(edgeData) {
    columns <- lapply(edgeData[.standardEdgeColumns()], as.character)
    do.call(paste, c(columns, sep = "\r"))
}

.fillBlank <- function(x, replacement) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- replacement
    x
}

.standardNodeColumns <- function() {
    c(
        "nodeKey",
        "nodeIdentifier",
        "nodeName",
        "assayName"
    )
}

.asPlainDataFrame <- function(x) {
    if (methods::is(x, "DataFrame")) {
        return(as.data.frame(x))
    }

    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

.coerceNumericMatrixStrict <- function(x, context = "assay matrix") {
    originalMatrix <- as.matrix(x)
    if (is.complex(originalMatrix) &&
        any(!is.na(originalMatrix) & Im(originalMatrix) != 0)) {
        stop("Complex values were found in ", context, ".", call. = FALSE)
    }
    originalMissing <- is.na(originalMatrix)
    numericMatrix <- originalMatrix
    suppressWarnings(storage.mode(numericMatrix) <- "numeric")
    dimnames(numericMatrix) <- dimnames(originalMatrix)

    introducedMissing <- is.na(numericMatrix) & !originalMissing
    if (any(introducedMissing)) {
        badIndex <- which(introducedMissing, arr.ind = TRUE)
        badIndex <- badIndex[seq_len(min(5L, nrow(badIndex))), , drop = FALSE]
        rowLabels <- rownames(originalMatrix)
        colLabels <- colnames(originalMatrix)
        badLocations <- vapply(
            seq_len(nrow(badIndex)),
            function(index) {
                rowLabel <- if (!is.null(rowLabels)) rowLabels[[badIndex[index, "row"]]] else badIndex[index, "row"]
                colLabel <- if (!is.null(colLabels)) colLabels[[badIndex[index, "col"]]] else badIndex[index, "col"]
                paste0(rowLabel, "/", colLabel)
            },
            character(1L)
        )
        stop(
            "Non-numeric values were found in ", context,
            ". First affected feature/sample locations: ",
            paste(badLocations, collapse = ", "),
            call. = FALSE
        )
    }

    nonFinite <- is.nan(numericMatrix) | is.infinite(numericMatrix)
    if (any(nonFinite)) {
        stop("Non-finite values were found in ", context, ".", call. = FALSE)
    }

    numericMatrix
}

.coerceNumericColumn <- function(x, columnName) {
    if (is.factor(x)) {
        x <- as.character(x)
    }
    if (is.complex(x) && any(!is.na(x) & Im(x) != 0)) {
        stop("Column `", columnName, "` cannot contain complex values.", call. = FALSE)
    }

    missingValue <- is.na(x)
    if (is.character(x)) {
        missingValue <- missingValue | !nzchar(trimws(x))
        x[missingValue] <- NA_character_
    }
    numericValue <- suppressWarnings(as.numeric(x))
    if (any(is.na(numericValue) & !missingValue)) {
        stop("Column `", columnName, "` must contain numeric values.", call. = FALSE)
    }
    if (any(is.nan(numericValue) | is.infinite(numericValue))) {
        stop("Column `", columnName, "` cannot contain non-finite values.", call. = FALSE)
    }

    numericValue
}

.coerceLogicalColumn <- function(x, columnName) {
    if (is.factor(x)) {
        x <- as.character(x)
    }
    missingValue <- is.na(x)
    if (is.character(x)) {
        missingValue <- missingValue | !nzchar(trimws(x))
        x[missingValue] <- NA_character_
    }
    logicalValue <- suppressWarnings(as.logical(x))
    if (any(is.na(logicalValue) & !missingValue)) {
        stop("Column `", columnName, "` must contain logical values.", call. = FALSE)
    }

    logicalValue
}

.hasAutomaticRowNames <- function(x) {
    rowNames <- rownames(x)
    if (is.null(rowNames)) {
        return(TRUE)
    }

    rowNamesInfo <- .row_names_info(x, type = 0L)
    is.integer(rowNamesInfo) &&
        length(rowNamesInfo) == 2L &&
        is.na(rowNamesInfo[1L]) &&
        rowNamesInfo[2L] < 0L
}

.emptyStandardNodeTable <- function(numberOfRows = 0L) {
    S4Vectors::DataFrame(
        nodeKey = rep(NA_character_, numberOfRows),
        nodeIdentifier = rep(NA_character_, numberOfRows),
        nodeName = rep(NA_character_, numberOfRows),
        assayName = rep(NA_character_, numberOfRows),
        check.names = FALSE
    )
}

.coerceStandardNodeTable <- function(nodeTable, preserveExtraColumns = FALSE) {
    if (is.null(nodeTable)) {
        return(.emptyStandardNodeTable())
    }

    nodeTable <- .asPlainDataFrame(nodeTable)
    nodeColumns <- .standardNodeColumns()
    extraColumns <- if (isTRUE(preserveExtraColumns)) {
        setdiff(names(nodeTable), nodeColumns)
    } else {
        character()
    }

    for (columnName in nodeColumns) {
        if (!columnName %in% names(nodeTable)) {
            nodeTable[[columnName]] <- rep(NA_character_, nrow(nodeTable))
        }
    }

    nodeTable <- nodeTable[, c(nodeColumns, extraColumns), drop = FALSE]
    for (columnName in nodeColumns) {
        nodeTable[[columnName]] <- as.character(nodeTable[[columnName]])
        nodeTable[[columnName]][is.na(nodeTable[[columnName]])] <- NA_character_
    }

    if (nrow(nodeTable)) {
        needsKey <- is.na(nodeTable$nodeKey) | !nzchar(nodeTable$nodeKey)
        nodeTable$nodeKey[needsKey] <- .nodeKey(
            assayName = nodeTable$assayName[needsKey],
            featureIdentifier = nodeTable$nodeIdentifier[needsKey]
        )

        needsName <- is.na(nodeTable$nodeName) | !nzchar(nodeTable$nodeName)
        nodeTable$nodeName[needsName] <- nodeTable$nodeIdentifier[needsName]

        nodeTable <- nodeTable[!duplicated(nodeTable$nodeKey), , drop = FALSE]
    }

    S4Vectors::DataFrame(nodeTable, check.names = FALSE)
}

.deriveNodeTableFromEdgeData <- function(edgeData) {
    if (!nrow(edgeData)) {
        return(.emptyStandardNodeTable())
    }

    nodeTable <- rbind(
        data.frame(
            nodeKey = .nodeKey(edgeData$fromAssayName, edgeData$fromFeatureIdentifier),
            nodeIdentifier = edgeData$fromFeatureIdentifier,
            nodeName = edgeData$fromFeatureName,
            assayName = edgeData$fromAssayName,
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        data.frame(
            nodeKey = .nodeKey(edgeData$toAssayName, edgeData$toFeatureIdentifier),
            nodeIdentifier = edgeData$toFeatureIdentifier,
            nodeName = edgeData$toFeatureName,
            assayName = edgeData$toAssayName,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    )

    .coerceStandardNodeTable(nodeTable)
}

.getStoredNetworkMetadata <- function(edgeTable) {
    if (!methods::is(edgeTable, "DataFrame")) {
        return(list())
    }

    storedMetadata <- S4Vectors::metadata(edgeTable)$cornetto
    if (!is.list(storedMetadata)) {
        return(list())
    }

    storedMetadata
}

.setStoredNetworkMetadata <- function(edgeTable, metadataList) {
    if (!methods::is(edgeTable, "DataFrame")) {
        edgeTable <- S4Vectors::DataFrame(edgeTable, check.names = FALSE)
    }

    storedMetadata <- .getStoredNetworkMetadata(edgeTable)
    if (length(metadataList)) {
        storedMetadata[names(metadataList)] <- metadataList
    }
    S4Vectors::metadata(edgeTable)$cornetto <- storedMetadata
    edgeTable
}

.getStoredNodeTable <- function(
    edgeTable,
    fallbackToEdges = FALSE,
    preserveExtraColumns = FALSE
) {
    storedMetadata <- .getStoredNetworkMetadata(edgeTable)
    if ("nodeTable" %in% names(storedMetadata)) {
        return(.coerceStandardNodeTable(
            storedMetadata$nodeTable,
            preserveExtraColumns = preserveExtraColumns
        ))
    }

    if (isTRUE(fallbackToEdges)) {
        edgeData <- .asPlainDataFrame(edgeTable)
        if (!all(c(
            "fromFeatureIdentifier",
            "toFeatureIdentifier",
            "fromFeatureName",
            "toFeatureName",
            "fromAssayName",
            "toAssayName"
        ) %in% names(edgeData))) {
            return(.emptyStandardNodeTable())
        }
        return(.deriveNodeTableFromEdgeData(edgeData))
    }

    NULL
}

.setStoredNodeTable <- function(edgeTable, nodeTable = NULL) {
    edgeTable <- if (methods::is(edgeTable, "DataFrame")) {
        edgeTable
    } else {
        S4Vectors::DataFrame(edgeTable, check.names = FALSE)
    }

    if (is.null(nodeTable)) {
        edgeData <- .asPlainDataFrame(edgeTable)
        nodeTable <- .deriveNodeTableFromEdgeData(edgeData)
    }

    .setStoredNetworkMetadata(
        edgeTable = edgeTable,
        metadataList = list(nodeTable = .coerceStandardNodeTable(nodeTable))
    )
}

.mergeNodeTables <- function(..., edgeTable = NULL) {
    inputs <- .flattenEdgeTables(...)
    nodeTables <- lapply(
        inputs,
        function(inputObject) {
            if (is.null(inputObject)) {
                return(NULL)
            }

            inputNames <- names(.asPlainDataFrame(inputObject))
            if (all(.standardNodeColumns() %in% inputNames)) {
                return(.coerceStandardNodeTable(inputObject))
            }

            .getStoredNodeTable(inputObject, fallbackToEdges = TRUE)
        }
    )

    if (!is.null(edgeTable)) {
        nodeTables <- c(nodeTables, list(.getStoredNodeTable(edgeTable, fallbackToEdges = TRUE)))
    }

    nodeTables <- Filter(
        function(nodeTable) {
            methods::is(nodeTable, "DataFrame") && nrow(nodeTable) > 0L
        },
        nodeTables
    )

    if (!length(nodeTables)) {
        return(.emptyStandardNodeTable())
    }

    .coerceStandardNodeTable(do.call(rbind, lapply(nodeTables, .asPlainDataFrame)))
}

.emptyStandardEdgeTable <- function(numberOfRows = 0L) {
    edgeColumns <- .standardEdgeColumns()
    numericColumns <- .numericEdgeColumns()
    logicalColumns <- .logicalEdgeColumns()
    tableList <- vector("list", length(edgeColumns))
    names(tableList) <- edgeColumns

    for (columnName in edgeColumns) {
        if (columnName %in% numericColumns) {
            tableList[[columnName]] <- rep(NA_real_, numberOfRows)
        } else if (columnName %in% logicalColumns) {
            tableList[[columnName]] <- rep(NA, numberOfRows)
        } else {
            tableList[[columnName]] <- rep(NA_character_, numberOfRows)
        }
    }

    .setStoredNodeTable(
        edgeTable = S4Vectors::DataFrame(tableList, check.names = FALSE),
        nodeTable = .emptyStandardNodeTable()
    )
}

.missingEdgeColumn <- function(columnName, numberOfRows) {
    if (columnName %in% .numericEdgeColumns()) {
        return(rep(NA_real_, numberOfRows))
    }
    if (columnName %in% .logicalEdgeColumns()) {
        return(rep(NA, numberOfRows))
    }

    rep(NA_character_, numberOfRows)
}

.coerceStandardEdgeTable <- function(edgeTable, nodeTable = NULL) {
    if (is.null(edgeTable)) {
        return(.emptyStandardEdgeTable())
    }

    preservedNodeTable <- nodeTable
    if (is.null(preservedNodeTable)) {
        preservedNodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = FALSE)
    }

    edgeTable <- .asPlainDataFrame(edgeTable)

    edgeColumns <- .standardEdgeColumns()
    numericColumns <- .numericEdgeColumns()
    logicalColumns <- .logicalEdgeColumns()

    for (columnName in edgeColumns) {
        if (!columnName %in% names(edgeTable)) {
            edgeTable[[columnName]] <- .missingEdgeColumn(
                columnName,
                nrow(edgeTable)
            )
        }
    }

    edgeTable <- edgeTable[, edgeColumns, drop = FALSE]

    for (columnName in numericColumns) {
        edgeTable[[columnName]] <- .coerceNumericColumn(
            edgeTable[[columnName]],
            columnName
        )
    }

    for (columnName in logicalColumns) {
        edgeTable[[columnName]] <- .coerceLogicalColumn(
            edgeTable[[columnName]],
            columnName
        )
    }

    remainingColumns <- setdiff(edgeColumns, c(numericColumns, logicalColumns))
    for (columnName in remainingColumns) {
        edgeTable[[columnName]] <- as.character(edgeTable[[columnName]])
        edgeTable[[columnName]][is.na(edgeTable[[columnName]])] <- NA_character_
    }

    .setStoredNodeTable(
        edgeTable = S4Vectors::DataFrame(edgeTable, check.names = FALSE),
        nodeTable = preservedNodeTable
    )
}

.matchDelimiter <- function(filePath, delimiter = NULL) {
    if (!is.null(delimiter)) {
        .assertScalarCharacter(delimiter, "delimiter")
        return(delimiter)
    }

    uncompressedPath <- sub(
        "[.](gz|bz2|xz)$",
        "",
        filePath,
        ignore.case = TRUE
    )
    extension <- tolower(tools::file_ext(uncompressedPath))
    if (identical(extension, "csv")) {
        return(",")
    }

    "\t"
}

.assertScalarCharacter <- function(x, argumentName) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop("`", argumentName, "` must be a single non-empty character value.", call. = FALSE)
    }
}

.assertScalarLogical <- function(x, argumentName) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("`", argumentName, "` must be TRUE or FALSE.", call. = FALSE)
    }
}

.validateUnitInterval <- function(x, argumentName, allowNull = FALSE) {
    if (is.null(x) && isTRUE(allowNull)) {
        return(NULL)
    }
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || x > 1) {
        stop(
            "`", argumentName, "` must be ",
            if (isTRUE(allowNull)) "NULL or " else "",
            "a numeric value between 0 and 1.",
            call. = FALSE
        )
    }

    as.numeric(x)
}

.validatePAdjustMethod <- function(pAdjustMethod) {
    .assertScalarCharacter(pAdjustMethod, "pAdjustMethod")
    allowed <- c(stats::p.adjust.methods, "qvalue")
    if (!pAdjustMethod %in% allowed) {
        stop(
            "`pAdjustMethod` must be one of: ",
            paste(allowed, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    pAdjustMethod
}

.assertWholeNumber <- function(x, argumentName, minimum = 0L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > .Machine$integer.max) {
        stop(
            "`", argumentName, "` must be a finite integer of at least ",
            minimum, ".",
            call. = FALSE
        )
    }

    as.integer(x)
}

# Fisher z is a large-sample approximation and its variance term is
# 1/(n1 - 3) + 1/(n2 - 3), so it degrades badly just above the n = 4 floor.
.warnOnSmallGroups <- function(groupLevels, group1Size, group2Size, threshold = 10L) {
    if (group1Size >= threshold && group2Size >= threshold) {
        return(invisible(NULL))
    }

    warning(
        "Differential correlation group sizes are below ", threshold,
        " for at least one group: ", groupLevels[[1L]], " = ", group1Size, ", ",
        groupLevels[[2L]], " = ", group2Size,
        ". The Fisher z difference is a large-sample approximation; below ",
        threshold, " samples per group the z-scores and p-values are unstable ",
        "and should be treated as descriptive.",
        call. = FALSE
    )
    invisible(NULL)
}

.assertTwoGroupLevels <- function(groupLevels) {
    groupLevels <- as.character(groupLevels)
    if (length(groupLevels) != 2L || anyNA(groupLevels) ||
        any(!nzchar(groupLevels)) || groupLevels[[1L]] == groupLevels[[2L]]) {
        stop(
            "`groupLevels` must contain two distinct, non-missing group labels.",
            call. = FALSE
        )
    }

    groupLevels
}

.assertMultiAssayExperiment <- function(analysisData) {
    if (!methods::is(analysisData, "MultiAssayExperiment")) {
        stop(
            "`analysisData` must be a MultiAssayExperiment object. ",
            "Use `createAnalysisData()` to build one from normalized assays.",
            call. = FALSE
        )
    }
}

.initializeCorNettoStore <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    store <- S4Vectors::metadata(analysisData)$cornetto
    if (is.null(store)) {
        store <- .corNettoStoreTemplate()
    } else {
        if (!is.list(store)) {
            stop(
                "`metadata(analysisData)$cornetto` is reserved for a named list of CorNetto results.",
                call. = FALSE
            )
        }
        template <- .corNettoStoreTemplate()
        missingSlots <- setdiff(names(template), names(store))
        if (length(missingSlots)) {
            store[missingSlots] <- template[missingSlots]
        }
    }
    analysisMetadata <- S4Vectors::metadata(analysisData)
    analysisMetadata[["cornetto"]] <- store
    S4Vectors::metadata(analysisData) <- analysisMetadata
    analysisData
}

.invalidateStoredResults <- function(analysisData) {
    store <- S4Vectors::metadata(analysisData)$cornetto
    if (is.list(store) && any(lengths(store) > 0L)) {
        warning(
            "Filtering invalidated stored CorNetto results; the result store was cleared.",
            call. = FALSE
        )
    }

    analysisMetadata <- S4Vectors::metadata(analysisData)
    analysisMetadata$cornetto <- .corNettoStoreTemplate()
    S4Vectors::metadata(analysisData) <- analysisMetadata
    analysisData
}

.storeCorNettoResults <- function(analysisData, slotName, resultObject, resultName) {
    analysisData <- .initializeCorNettoStore(analysisData)
    .assertScalarCharacter(slotName, "slotName")
    .assertScalarCharacter(resultName, "resultName")

    store <- S4Vectors::metadata(analysisData)$cornetto
    if (!slotName %in% names(store)) {
        stop("Unknown CorNetto result slot: `", slotName, "`.", call. = FALSE)
    }

    resultSlot <- store[[slotName]]
    resultSlot[[resultName]] <- resultObject
    store[[slotName]] <- resultSlot
    analysisMetadata <- S4Vectors::metadata(analysisData)
    analysisMetadata[["cornetto"]] <- store
    S4Vectors::metadata(analysisData) <- analysisMetadata
    analysisData
}

.makeResultName <- function(..., separator = "__") {
    pieces <- unlist(list(...), use.names = FALSE)
    pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
    paste(pieces, collapse = separator)
}

## assay(se) with no `i` silently returns the first matrix, which is the wrong
## answer whenever a user hands us counts alongside normalized values.
.selectAssayMatrix <- function(assayObject, assayName) {
    availableAssays <- SummarizedExperiment::assayNames(assayObject)
    if ("abundance" %in% availableAssays) {
        return(SummarizedExperiment::assay(assayObject, "abundance"))
    }
    if (length(availableAssays) == 1L) {
        return(SummarizedExperiment::assay(assayObject, availableAssays[[1L]]))
    }

    stop(
        "Assay `", assayName, "` is a SummarizedExperiment holding ",
        length(availableAssays), " matrices (",
        paste(availableAssays, collapse = ", "),
        "). CorNetto correlates one normalized matrix per assay: name it ",
        "`abundance`, or subset the object before calling createAnalysisData().",
        call. = FALSE
    )
}

.extractAssayMatrix <- function(analysisData, assayName) {
    .assertMultiAssayExperiment(analysisData)
    .assertScalarCharacter(assayName, "assayName")
    availableAssays <- names(MultiAssayExperiment::experiments(analysisData))
    if (!assayName %in% availableAssays) {
        stop(
            "`assayName` was not found in `analysisData`. Available assays are: ",
            paste(availableAssays, collapse = ", "),
            call. = FALSE
        )
    }

    assayObject <- MultiAssayExperiment::experiments(analysisData)[[assayName]]
    if (methods::is(assayObject, "SummarizedExperiment")) {
        assayMatrix <- .selectAssayMatrix(assayObject, assayName)
    } else if (is.matrix(assayObject) || is.data.frame(assayObject)) {
        assayMatrix <- as.matrix(assayObject)
    } else {
        stop(
            "Unsupported assay class for `", assayName, "`: ",
            paste(class(assayObject), collapse = ", "),
            call. = FALSE
        )
    }

    .coerceNumericMatrixStrict(
        assayMatrix,
        context = paste0("assay `", assayName, "`")
    )
}

.extractFeatureAnnotations <- function(analysisData, assayName, featureNameColumn = "featureName") {
    assayObject <- MultiAssayExperiment::experiments(analysisData)[[assayName]]
    featureIdentifiers <- rownames(.extractAssayMatrix(analysisData, assayName))
    featureNames <- featureIdentifiers

    if (methods::is(assayObject, "SummarizedExperiment")) {
        rowDataFrame <- .asPlainDataFrame(SummarizedExperiment::rowData(assayObject))
        if (featureNameColumn %in% names(rowDataFrame)) {
            featureNames <- as.character(rowDataFrame[[featureNameColumn]])
            featureNames[is.na(featureNames) | !nzchar(featureNames)] <- featureIdentifiers[
                is.na(featureNames) | !nzchar(featureNames)
            ]
        }
    }

    S4Vectors::DataFrame(
        featureIdentifier = featureIdentifiers,
        featureName = featureNames,
        assayName = assayName,
        check.names = FALSE
    )
}

.resolveFeatureSubset <- function(analysisData, assayName, featureSubset = NULL) {
    assayMatrix <- .extractAssayMatrix(analysisData, assayName)
    if (is.null(featureSubset)) {
        return(rownames(assayMatrix))
    }

    featureSubset <- unique(as.character(featureSubset))
    if (anyNA(featureSubset) || any(!nzchar(featureSubset))) {
        stop("`featureSubset` must be non-missing and non-empty.", call. = FALSE)
    }
    missingFeatures <- setdiff(featureSubset, rownames(assayMatrix))
    if (length(missingFeatures)) {
        stop(
            "Requested features were not found in assay `", assayName, "`: ",
            paste(utils::head(missingFeatures, 5L), collapse = ", "),
            call. = FALSE
        )
    }

    featureSubset
}

.resolveGroupSamples <- function(
    analysisData,
    groupColumn = NULL,
    groupLevel = NULL,
    sampleIds = NULL
) {
    sampleData <- .asPlainDataFrame(S4Vectors::DataFrame(MultiAssayExperiment::colData(analysisData)))
    if (is.null(rownames(sampleData))) {
        stop("`analysisData` must have row names in `colData`.", call. = FALSE)
    }

    if (!is.null(sampleIds)) {
        if (!is.null(groupColumn) || !is.null(groupLevel)) {
            stop(
                "Supply either `sampleIds` or `groupColumn` and `groupLevel`, not both.",
                call. = FALSE
            )
        }
        sampleIds <- as.character(sampleIds)
        if (anyNA(sampleIds) || any(!nzchar(sampleIds))) {
            stop("`sampleIds` must be non-missing and non-empty.", call. = FALSE)
        }
        missingSampleIds <- setdiff(sampleIds, rownames(sampleData))
        if (length(missingSampleIds)) {
            stop(
                "`sampleIds` were not found in `analysisData`: ",
                paste(utils::head(missingSampleIds, 5L), collapse = ", "),
                call. = FALSE
            )
        }
        return(intersect(rownames(sampleData), sampleIds))
    }

    if (is.null(groupColumn) && is.null(groupLevel)) {
        return(rownames(sampleData))
    }
    if (is.null(groupColumn) || is.null(groupLevel)) {
        stop(
            "`groupColumn` and `groupLevel` must be supplied together ",
            "unless `sampleIds` is supplied.",
            call. = FALSE
        )
    }

    .assertScalarCharacter(groupColumn, "groupColumn")
    .assertScalarCharacter(groupLevel, "groupLevel")

    if (!groupColumn %in% names(sampleData)) {
        stop(
            "`groupColumn` was not found in `colData`: ",
            groupColumn,
            call. = FALSE
        )
    }

    sampleIds <- rownames(sampleData)[sampleData[[groupColumn]] %in% groupLevel]
    if (!length(sampleIds)) {
        stop(
            "No samples were found for `groupLevel = ", groupLevel,
            "` in `groupColumn = ", groupColumn, "`.",
            call. = FALSE
        )
    }

    sampleIds
}

.flattenEdgeTables <- function(...) {
    inputs <- list(...)
    flattened <- list()

    for (inputObject in inputs) {
        if (is.null(inputObject)) {
            next
        }

        if (is.list(inputObject) &&
            !is.data.frame(inputObject) &&
            !methods::is(inputObject, "DataFrame")) {
            flattened <- c(flattened, unname(inputObject))
        } else {
            flattened <- c(flattened, list(inputObject))
        }
    }

    flattened
}

.coerceSampleData <- function(sampleData) {
    if (is.null(sampleData)) {
        stop("`sampleData` must not be NULL.", call. = FALSE)
    }

    sampleData <- .asPlainDataFrame(sampleData)

    if (.hasAutomaticRowNames(sampleData)) {
        if (!"sampleId" %in% names(sampleData)) {
            stop(
                "`sampleData` must have explicit row names or a `sampleId` column.",
                call. = FALSE
            )
        }
        rownames(sampleData) <- as.character(sampleData$sampleId)
    }

    if (anyNA(rownames(sampleData)) || any(!nzchar(rownames(sampleData))) ||
        anyDuplicated(rownames(sampleData))) {
        stop("`sampleData` row names must be unique, non-missing sample identifiers.", call. = FALSE)
    }

    S4Vectors::DataFrame(sampleData, check.names = FALSE)
}

#' Get Stored CorNetto Results
#'
#' Retrieve the complete CorNetto result store, a result-type slot, or a
#' single named stored result from a `MultiAssayExperiment`.
#'
#' @param analysisData A `MultiAssayExperiment` object.
#' @param resultType Optional result type. When `NULL`, the complete
#'   result store is returned unless `resultName` is supplied. Common
#'   values include `"knowledgeNetworks"`, `"correlationResults"`,
#'   `"differentialCorrelationResults"`,
#'   `"differentialCorrelationNetworks"`, `"rewiringResults"`,
#'   `"validationResults"`, and `"integratedNetworks"`.
#' @param resultName Optional stored result name. When supplied with
#'   `resultType`, a single result is returned. When supplied without
#'   `resultType`, all result slots are searched and the result name must
#'   be unique.
#'
#' @return A named list of stored CorNetto results, a single result-type
#'   list, or one stored result object.
#' @examples
#' analysisData <- exampleAnalysisData()
#' analysisData <- validateKnowledgeNetwork(
#'     exampleKnowledgeNetwork(),
#'     analysisData = analysisData,
#'     storeResult = TRUE
#' )
#' names(getCorNettoResult(analysisData))
#' names(knowledgeNetworks(analysisData))
#' getCorNettoResult(analysisData, "knowledgeNetworks", "knowledgeNetwork")
#' @export
getCorNettoResult <- function(
    analysisData,
    resultType = NULL,
    resultName = NULL
) {
    .assertMultiAssayExperiment(analysisData)
    if (!is.null(resultType)) {
        .assertScalarCharacter(resultType, "resultType")
    }
    if (!is.null(resultName)) {
        .assertScalarCharacter(resultName, "resultName")
    }
    analysisData <- .initializeCorNettoStore(analysisData)
    store <- S4Vectors::metadata(analysisData)$cornetto

    if (is.null(resultType)) {
        if (is.null(resultName)) {
            return(store)
        }
        matchingTypes <- names(Filter(
            function(resultSlot) resultName %in% names(resultSlot),
            store
        ))
        if (!length(matchingTypes)) {
            stop("No stored result named `", resultName, "` was found.", call. = FALSE)
        }
        if (length(matchingTypes) > 1L) {
            stop(
                "`resultName` was found in multiple result types: ",
                paste(matchingTypes, collapse = ", "),
                ". Supply `resultType` to disambiguate.",
                call. = FALSE
            )
        }
        return(store[[matchingTypes[[1L]]]][[resultName]])
    }

    if (!resultType %in% names(store)) {
        stop("Unknown `resultType`: ", resultType, call. = FALSE)
    }

    resultSlot <- store[[resultType]]
    if (is.null(resultName)) {
        return(resultSlot)
    }
    if (!resultName %in% names(resultSlot)) {
        stop(
            "No stored result named `", resultName,
            "` was found in result type `", resultType, "`.",
            call. = FALSE
        )
    }
    resultSlot[[resultName]]
}

#' @rdname getCorNettoResult
#' @export
knowledgeNetworks <- function(analysisData) {
    getCorNettoResult(analysisData, "knowledgeNetworks")
}

#' @rdname getCorNettoResult
#' @export
correlationResults <- function(analysisData) {
    getCorNettoResult(analysisData, "correlationResults")
}

#' @rdname getCorNettoResult
#' @export
differentialCorrelationResults <- function(analysisData) {
    getCorNettoResult(analysisData, "differentialCorrelationResults")
}

#' @rdname getCorNettoResult
#' @export
differentialCorrelationNetworks <- function(analysisData) {
    getCorNettoResult(analysisData, "differentialCorrelationNetworks")
}

#' @rdname getCorNettoResult
#' @export
rewiringResults <- function(analysisData) {
    getCorNettoResult(analysisData, "rewiringResults")
}

#' @rdname getCorNettoResult
#' @export
validationResults <- function(analysisData) {
    getCorNettoResult(analysisData, "validationResults")
}

#' @rdname getCorNettoResult
#' @export
integratedNetworks <- function(analysisData) {
    getCorNettoResult(analysisData, "integratedNetworks")
}
