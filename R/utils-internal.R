.corNettoStoreTemplate <- function() {
    list(
        knowledgeNetworks = list(),
        correlationNetworks = list(),
        sparseWithinOmicCorrelations = list(),
        sparseCrossOmicCorrelations = list(),
        differentialCorrelationResults = list(),
        rewiringResults = list(),
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
        "edgeWeight",
        "evidenceScore"
    )
}

.logicalEdgeColumns <- function() {
    "isDirected"
}

.standardNodeColumns <- function() {
    c(
        "nodeKey",
        "nodeIdentifier",
        "nodeName",
        "assayName"
    )
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

    nodeTable <- as.data.frame(nodeTable, stringsAsFactors = FALSE, check.names = FALSE)
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
        edgeData <- as.data.frame(edgeTable, stringsAsFactors = FALSE, check.names = FALSE)
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
        edgeData <- as.data.frame(edgeTable, stringsAsFactors = FALSE, check.names = FALSE)
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

            inputNames <- names(as.data.frame(inputObject, stringsAsFactors = FALSE, check.names = FALSE))
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

    .coerceStandardNodeTable(do.call(rbind, lapply(nodeTables, as.data.frame, stringsAsFactors = FALSE)))
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
        edgeTable <- as.data.frame(edgeTable, stringsAsFactors = FALSE)
    } else {
        edgeTable <- as.data.frame(edgeTable, stringsAsFactors = FALSE)
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
        S4Vectors::metadata(analysisData)$cornetto <- .corNettoStoreTemplate()
    }
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

    store[[slotName]][[resultName]] <- resultObject
    S4Vectors::metadata(analysisData)$cornetto <- store
    analysisData
}

.getCorNettoResults <- function(analysisData, slotName) {
    analysisData <- .initializeCorNettoStore(analysisData)
    store <- S4Vectors::metadata(analysisData)$cornetto
    store[[slotName]]
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

    storage.mode(assayMatrix) <- "numeric"
    assayMatrix
}

.extractFeatureAnnotations <- function(analysisData, assayName, featureNameColumn = "featureName") {
    assayObject <- MultiAssayExperiment::experiments(analysisData)[[assayName]]
    featureIdentifiers <- rownames(.extractAssayMatrix(analysisData, assayName))
    featureNames <- featureIdentifiers

    if (methods::is(assayObject, "SummarizedExperiment")) {
        rowDataFrame <- as.data.frame(SummarizedExperiment::rowData(assayObject))
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
    sampleData <- as.data.frame(S4Vectors::DataFrame(MultiAssayExperiment::colData(analysisData)))
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
        sampleData <- as.data.frame(sampleData, stringsAsFactors = FALSE)
    } else {
        sampleData <- as.data.frame(sampleData, stringsAsFactors = FALSE)
    }

    if (is.null(rownames(sampleData))) {
        if ("sampleId" %in% names(sampleData)) {
            rownames(sampleData) <- as.character(sampleData$sampleId)
        } else {
            stop(
                "`sampleData` must have row names or a `sampleId` column.",
                call. = FALSE
            )
        }
    }

    S4Vectors::DataFrame(sampleData, check.names = FALSE)
}

.isObservedAssayPair <- function(edgeTable) {
    edgeTable$fromAssayName %in% edgeTable$toAssayName
}

#' CorNetto Results
#'
#' Retrieve the complete CorNetto result store attached to a
#' `MultiAssayExperiment`.
#'
#' @param analysisData A `MultiAssayExperiment` object.
#'
#' @return A named list of stored CorNetto results.
#' @export
corNettoResults <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    analysisData <- .initializeCorNettoStore(analysisData)
    S4Vectors::metadata(analysisData)$cornetto
}

#' @export
knowledgeNetworks <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "knowledgeNetworks")
}

#' @export
correlationNetworks <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "correlationNetworks")
}

#' @export
sparseWithinOmicCorrelations <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "sparseWithinOmicCorrelations")
}

#' @export
sparseCrossOmicCorrelations <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "sparseCrossOmicCorrelations")
}

#' @export
differentialCorrelationResults <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "differentialCorrelationResults")
}

#' @export
rewiringResults <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "rewiringResults")
}

#' @export
integratedNetworks <- function(analysisData) {
    .assertMultiAssayExperiment(analysisData)
    .getCorNettoResults(.initializeCorNettoStore(analysisData), "integratedNetworks")
}
