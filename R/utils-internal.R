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

    nonFinite <- !is.na(numericMatrix) & !is.finite(numericMatrix)
    if (any(nonFinite)) {
        stop("Non-finite values were found in ", context, ".", call. = FALSE)
    }

    numericMatrix
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

.coerceStandardNodeTable <- function(nodeTable) {
    if (is.null(nodeTable)) {
        return(.emptyStandardNodeTable())
    }

    nodeTable <- .asPlainDataFrame(nodeTable)
    nodeColumns <- .standardNodeColumns()

    for (columnName in nodeColumns) {
        if (!columnName %in% names(nodeTable)) {
            nodeTable[[columnName]] <- NA_character_
        }
    }

    nodeTable <- nodeTable[, nodeColumns, drop = FALSE]
    for (columnName in nodeColumns) {
        nodeTable[[columnName]] <- as.character(nodeTable[[columnName]])
        nodeTable[[columnName]][is.na(nodeTable[[columnName]])] <- NA_character_
    }

    if (nrow(nodeTable)) {
        nodeTable$nodeKey[
            is.na(nodeTable$nodeKey) | !nzchar(nodeTable$nodeKey)
        ] <- .nodeKey(
            assayName = nodeTable$assayName[
                is.na(nodeTable$nodeKey) | !nzchar(nodeTable$nodeKey)
            ],
            featureIdentifier = nodeTable$nodeIdentifier[
                is.na(nodeTable$nodeKey) | !nzchar(nodeTable$nodeKey)
            ]
        )
        nodeTable$nodeName[
            is.na(nodeTable$nodeName) | !nzchar(nodeTable$nodeName)
        ] <- nodeTable$nodeIdentifier[
            is.na(nodeTable$nodeName) | !nzchar(nodeTable$nodeName)
        ]
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

.getStoredNodeTable <- function(edgeTable, fallbackToEdges = FALSE) {
    storedMetadata <- .getStoredNetworkMetadata(edgeTable)
    if ("nodeTable" %in% names(storedMetadata)) {
        return(.coerceStandardNodeTable(storedMetadata$nodeTable))
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

.coerceStandardEdgeTable <- function(edgeTable, nodeTable = NULL) {
    if (is.null(edgeTable)) {
        return(.emptyStandardEdgeTable())
    }

    preservedNodeTable <- nodeTable
    if (is.null(preservedNodeTable)) {
        preservedNodeTable <- .getStoredNodeTable(edgeTable, fallbackToEdges = FALSE)
    }

    if (methods::is(edgeTable, "DataFrame")) {
        edgeTable <- .asPlainDataFrame(edgeTable)
    } else {
        edgeTable <- .asPlainDataFrame(edgeTable)
    }

    edgeColumns <- .standardEdgeColumns()
    numericColumns <- .numericEdgeColumns()
    logicalColumns <- .logicalEdgeColumns()

    for (columnName in edgeColumns) {
        if (!columnName %in% names(edgeTable)) {
            if (columnName %in% numericColumns) {
                edgeTable[[columnName]] <- NA_real_
            } else if (columnName %in% logicalColumns) {
                edgeTable[[columnName]] <- NA
            } else {
                edgeTable[[columnName]] <- NA_character_
            }
        }
    }

    edgeTable <- edgeTable[, edgeColumns, drop = FALSE]

    for (columnName in numericColumns) {
        edgeTable[[columnName]] <- as.numeric(edgeTable[[columnName]])
    }

    for (columnName in logicalColumns) {
        edgeTable[[columnName]] <- as.logical(edgeTable[[columnName]])
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
        return(delimiter)
    }

    extension <- tolower(tools::file_ext(filePath))
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

.getCorNettoResults <- function(analysisData, slotName) {
    analysisData <- .initializeCorNettoStore(analysisData)
    store <- S4Vectors::metadata(analysisData)$cornetto
    store[[slotName]]
}

.resolveResultType <- function(resultType) {
    resultAliases <- c(
        combinedNetworks = "integratedNetworks"
    )
    if (resultType %in% names(resultAliases)) {
        return(unname(resultAliases[[resultType]]))
    }
    resultType
}

.makeResultName <- function(..., separator = "__") {
    pieces <- unlist(list(...), use.names = FALSE)
    pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
    paste(pieces, collapse = separator)
}

.listToNamedVector <- function(x, expectedLength = NULL) {
    if (is.null(x)) {
        return(NULL)
    }

    if (is.list(x) && !methods::is(x, "DataFrame")) {
        return(x)
    }

    if (!is.null(expectedLength)) {
        rep(list(x), expectedLength)
    } else {
        list(x)
    }
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
        assayMatrix <- SummarizedExperiment::assay(assayObject)
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

    featureSubset <- as.character(featureSubset)
    featureSubset <- intersect(featureSubset, rownames(assayMatrix))
    if (!length(featureSubset)) {
        stop(
            "No requested features were found in assay `", assayName, "`.",
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
        sampleIds <- as.character(sampleIds)
        sampleIds <- intersect(sampleIds, rownames(sampleData))
        if (!length(sampleIds)) {
            stop("None of the requested `sampleIds` were found in `analysisData`.", call. = FALSE)
        }
        return(sampleIds)
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

        if (is.list(inputObject) && !methods::is(inputObject, "DataFrame")) {
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

    if (methods::is(sampleData, "DataFrame")) {
        sampleData <- .asPlainDataFrame(sampleData)
    } else {
        sampleData <- .asPlainDataFrame(sampleData)
    }

    if ("sampleId" %in% names(sampleData) && .hasAutomaticRowNames(sampleData)) {
        rownames(sampleData) <- as.character(sampleData$sampleId)
    }

    if (is.null(rownames(sampleData))) {
        stop(
            "`sampleData` must have row names or a `sampleId` column.",
            call. = FALSE
        )
    }

    if (anyNA(rownames(sampleData)) || any(!nzchar(rownames(sampleData))) ||
        anyDuplicated(rownames(sampleData))) {
        stop("`sampleData` row names must be unique, non-missing sample identifiers.", call. = FALSE)
    }

    S4Vectors::DataFrame(sampleData, check.names = FALSE)
}

.isObservedAssayPair <- function(edgeTable) {
    edgeTable$fromAssayName == edgeTable$toAssayName
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

    resultType <- .resolveResultType(resultType)
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
