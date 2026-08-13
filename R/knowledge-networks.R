.standardizeKnowledgeNetwork <- function(
    knowledgeNetwork,
    edgeType = "association",
    edgeDirection = "undirected",
    knowledgeSource = "userSupplied"
) {
    knowledgeNetwork <- .asPlainDataFrame(knowledgeNetwork)
    requiredColumns <- c(
        "fromFeatureIdentifier",
        "toFeatureIdentifier",
        "fromAssayName",
        "toAssayName"
    )

    missingColumns <- setdiff(requiredColumns, names(knowledgeNetwork))
    if (length(missingColumns)) {
        stop(
            "Knowledge network is missing required columns: ",
            paste(missingColumns, collapse = ", "),
            call. = FALSE
        )
    }

    if (!"edgeType" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$edgeType <- edgeType
    }
    if (!"edgeDirection" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$edgeDirection <- edgeDirection
    }
    if (!"knowledgeSource" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$knowledgeSource <- knowledgeSource
    }
    if (!"evidenceScore" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$evidenceScore <- NA_real_
    }
    knowledgeNetwork$evidenceScore <- .coerceNumericColumn(
        knowledgeNetwork$evidenceScore,
        "evidenceScore"
    )
    if (!"fromFeatureName" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$fromFeatureName <- knowledgeNetwork$fromFeatureIdentifier
    }
    if (!"toFeatureName" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$toFeatureName <- knowledgeNetwork$toFeatureIdentifier
    }

    # A blank cell in a column that exists must fall back to the same default as
    # a missing column, otherwise isDirected() below reads it as "not undirected"
    # and silently makes the edge directed.
    knowledgeNetwork$edgeType <- .fillBlank(knowledgeNetwork$edgeType, edgeType)
    knowledgeNetwork$edgeDirection <- .fillBlank(knowledgeNetwork$edgeDirection, edgeDirection)
    knowledgeNetwork$knowledgeSource <- .fillBlank(knowledgeNetwork$knowledgeSource, knowledgeSource)

    if (!"isDirected" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$isDirected <-
            !(tolower(knowledgeNetwork$edgeDirection) %in% c("undirected", "association")) &
            !(tolower(knowledgeNetwork$edgeType) %in% c("proteinproteininteraction", "correlation"))
    } else {
        # An explicit column wins, but a blank entry falls back to inference.
        suppliedDirection <- as.logical(knowledgeNetwork$isDirected)
        inferredDirection <-
            !(tolower(knowledgeNetwork$edgeDirection) %in% c("undirected", "association")) &
            !(tolower(knowledgeNetwork$edgeType) %in% c("proteinproteininteraction", "correlation"))
        suppliedDirection[is.na(suppliedDirection)] <- inferredDirection[is.na(suppliedDirection)]
        knowledgeNetwork$isDirected <- suppliedDirection
    }

    knowledgeNetwork$sourceType <- "knowledge"
    knowledgeNetwork$correlationScope <- NA_character_
    knowledgeNetwork$correlationMethod <- NA_character_
    knowledgeNetwork$groupName <- NA_character_
    knowledgeNetwork$comparisonName <- NA_character_
    knowledgeNetwork$correlationValue <- NA_real_
    knowledgeNetwork$group1CorrelationValue <- NA_real_
    knowledgeNetwork$group2CorrelationValue <- NA_real_
    knowledgeNetwork$pValue <- NA_real_
    knowledgeNetwork$adjustedPValue <- NA_real_
    knowledgeNetwork$group1PValue <- NA_real_
    knowledgeNetwork$group2PValue <- NA_real_
    knowledgeNetwork$group1AdjustedPValue <- NA_real_
    knowledgeNetwork$group2AdjustedPValue <- NA_real_
    knowledgeNetwork$zScoreDifference <- NA_real_
    knowledgeNetwork$edgeWeight <- ifelse(
        is.na(knowledgeNetwork$evidenceScore),
        1,
        knowledgeNetwork$evidenceScore
    )

    .coerceStandardEdgeTable(knowledgeNetwork)
}

#' Validate a Knowledge Network
#'
#' Validate and standardize a user-supplied prior-knowledge network.
#'
#' @param knowledgeNetwork A data frame-like object containing prior
#'   interactions.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   standardized network.
#' @param resultName Name used when storing the knowledge network in
#'   `metadata(analysisData)$cornetto`.
#' @param storeResult Whether to store the result in `analysisData`.
#' @param edgeType Default edge type when the column is missing.
#' @param edgeDirection Default edge direction when the column is
#'   missing.
#' @param knowledgeSource Default source label when the column is
#'   missing.
#' @param quiet Whether to suppress the success message.
#'
#' @return A standardized `DataFrame` or an updated
#'   `MultiAssayExperiment` when `storeResult = TRUE`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' validated <- validateKnowledgeNetwork(knowledgeNetwork)
#' validated
#' @export
validateKnowledgeNetwork <- function(
    knowledgeNetwork,
    analysisData = NULL,
    resultName = "knowledgeNetwork",
    storeResult = FALSE,
    edgeType = "association",
    edgeDirection = "undirected",
    knowledgeSource = "userSupplied",
    quiet = FALSE
) {
    .assertScalarLogical(quiet, "quiet")
    .assertScalarLogical(storeResult, "storeResult")
    .assertScalarCharacter(edgeType, "edgeType")
    .assertScalarCharacter(edgeDirection, "edgeDirection")
    .assertScalarCharacter(knowledgeSource, "knowledgeSource")
    standardizedNetwork <- .standardizeKnowledgeNetwork(
        knowledgeNetwork = knowledgeNetwork,
        edgeType = edgeType,
        edgeDirection = edgeDirection,
        knowledgeSource = knowledgeSource
    )

    if (anyNA(standardizedNetwork$fromFeatureIdentifier) ||
        anyNA(standardizedNetwork$toFeatureIdentifier) ||
        any(!nzchar(trimws(standardizedNetwork$fromFeatureIdentifier))) ||
        any(!nzchar(trimws(standardizedNetwork$toFeatureIdentifier)))) {
        stop(
            "Knowledge networks cannot contain missing or blank feature identifiers.",
            call. = FALSE
        )
    }
    if (anyNA(standardizedNetwork$fromAssayName) ||
        anyNA(standardizedNetwork$toAssayName) ||
        any(!nzchar(trimws(standardizedNetwork$fromAssayName))) ||
        any(!nzchar(trimws(standardizedNetwork$toAssayName)))) {
        stop(
            "Knowledge networks cannot contain missing or blank assay names.",
            call. = FALSE
        )
    }

    if (!storeResult) {
        if (!isTRUE(quiet)) {
            message(
                "Knowledge network is structured correctly for CorNetto."
            )
        }
        return(standardizedNetwork)
    }

    .assertMultiAssayExperiment(analysisData)
    analysisData <- .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "knowledgeNetworks",
        resultObject = standardizedNetwork,
        resultName = resultName
    )
    if (!isTRUE(quiet)) {
        message("Knowledge network is structured correctly for CorNetto.")
    }
    analysisData
}

#' Summarize a Knowledge Network
#'
#' Summarize prior-network size, edge annotations, assay-pair coverage,
#' duplicate edges, and optionally how much of the network is measurable
#' in a CorNetto analysis object.
#'
#' @param knowledgeNetwork A standardized knowledge network or any table
#'   accepted by `validateKnowledgeNetwork()`.
#' @param analysisData Optional `MultiAssayExperiment` used to count
#'   edges whose assay and feature endpoints are present in measured
#'   data.
#' @param quiet Whether to suppress emitted warnings and only return
#'   them in the result object.
#'
#' @return A named list with summary `DataFrame` objects and warnings.
#' @examples
#' analysisData <- exampleAnalysisData()
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' summary <- summarizeKnowledgeNetwork(knowledgeNetwork, analysisData)
#' summary$overall
#' summary$assayPairs
#' @export
summarizeKnowledgeNetwork <- function(
    knowledgeNetwork,
    analysisData = NULL,
    quiet = FALSE
) {
    .assertScalarLogical(quiet, "quiet")
    knowledgeNetwork <- validateKnowledgeNetwork(
        knowledgeNetwork,
        storeResult = FALSE,
        quiet = TRUE
    )
    edgeData <- .asPlainDataFrame(knowledgeNetwork)
    nodeTable <- .getStoredNodeTable(knowledgeNetwork, fallbackToEdges = TRUE)
    deduplicatedEdges <- .deduplicateUndirectedEdges(knowledgeNetwork)
    warningMessages <- character()

    if (!nrow(edgeData)) {
        warningMessages <- c(warningMessages, "Knowledge network contains no edges.")
    }

    duplicateEdgeCount <- nrow(edgeData) - nrow(deduplicatedEdges)
    if (duplicateEdgeCount > 0L) {
        warningMessages <- c(
            warningMessages,
            paste0("Knowledge network contains ", duplicateEdgeCount, " duplicated edges.")
        )
    }

    measuredAssayEdgeCount <- NA_integer_
    measuredFeatureEdgeCount <- NA_integer_
    unmeasuredAssayEdgeCount <- NA_integer_
    unmeasuredFeatureEdgeCount <- NA_integer_

    if (!is.null(analysisData)) {
        analysisData <- validateAnalysisData(analysisData, quiet = TRUE)
        assayNames <- names(MultiAssayExperiment::experiments(analysisData))
        featureSets <- lapply(
            assayNames,
            function(assayName) rownames(.extractAssayMatrix(analysisData, assayName))
        )
        names(featureSets) <- assayNames

        fromAssayMeasured <- edgeData$fromAssayName %in% assayNames
        toAssayMeasured <- edgeData$toAssayName %in% assayNames
        measuredAssayEdge <- fromAssayMeasured & toAssayMeasured

        fromFeatureMeasured <- mapply(
            function(assayName, featureIdentifier) {
                assayName %in% names(featureSets) &&
                    featureIdentifier %in% featureSets[[assayName]]
            },
            edgeData$fromAssayName,
            edgeData$fromFeatureIdentifier
        )
        toFeatureMeasured <- mapply(
            function(assayName, featureIdentifier) {
                assayName %in% names(featureSets) &&
                    featureIdentifier %in% featureSets[[assayName]]
            },
            edgeData$toAssayName,
            edgeData$toFeatureIdentifier
        )
        measuredFeatureEdge <- measuredAssayEdge & fromFeatureMeasured & toFeatureMeasured

        measuredAssayEdgeCount <- sum(measuredAssayEdge)
        measuredFeatureEdgeCount <- sum(measuredFeatureEdge)
        unmeasuredAssayEdgeCount <- nrow(edgeData) - measuredAssayEdgeCount
        unmeasuredFeatureEdgeCount <- nrow(edgeData) - measuredFeatureEdgeCount

        if (measuredFeatureEdgeCount == 0L && nrow(edgeData) > 0L) {
            warningMessages <- c(
                warningMessages,
                "No knowledge-network edges have both endpoints measured in `analysisData`."
            )
        }
    }

    missingFeatureNameCount <- sum(
        is.na(edgeData$fromFeatureName) | !nzchar(edgeData$fromFeatureName) |
            is.na(edgeData$toFeatureName) | !nzchar(edgeData$toFeatureName)
    )

    evidenceValues <- as.numeric(edgeData$evidenceScore)
    evidenceValues <- evidenceValues[!is.na(evidenceValues)]
    evidenceSummary <- if (length(evidenceValues)) {
        data.frame(
            nonMissingEvidenceScores = length(evidenceValues),
            minimumEvidenceScore = min(evidenceValues),
            medianEvidenceScore = stats::median(evidenceValues),
            maximumEvidenceScore = max(evidenceValues),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    } else {
        data.frame(
            nonMissingEvidenceScores = 0L,
            minimumEvidenceScore = NA_real_,
            medianEvidenceScore = NA_real_,
            maximumEvidenceScore = NA_real_,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    assayPair <- paste(edgeData$fromAssayName, edgeData$toAssayName, sep = " -> ")
    assayPairTable <- as.data.frame(table(assayPair), stringsAsFactors = FALSE)
    names(assayPairTable) <- c("assayPair", "edgeCount")

    edgeTypeTable <- as.data.frame(table(edgeData$edgeType), stringsAsFactors = FALSE)
    names(edgeTypeTable) <- c("edgeType", "edgeCount")

    edgeDirectionTable <- as.data.frame(table(edgeData$edgeDirection), stringsAsFactors = FALSE)
    names(edgeDirectionTable) <- c("edgeDirection", "edgeCount")

    sourceTable <- as.data.frame(table(edgeData$knowledgeSource), stringsAsFactors = FALSE)
    names(sourceTable) <- c("knowledgeSource", "edgeCount")

    warningMessages <- unique(warningMessages)
    if (!isTRUE(quiet) && length(warningMessages)) {
        warning(paste(warningMessages, collapse = "\n"), call. = FALSE)
    }

    list(
        overall = S4Vectors::DataFrame(
            totalEdges = nrow(edgeData),
            uniqueNodes = nrow(nodeTable),
            duplicateEdges = duplicateEdgeCount,
            missingFeatureNameEdges = missingFeatureNameCount,
            measuredAssayEdges = measuredAssayEdgeCount,
            unmeasuredAssayEdges = unmeasuredAssayEdgeCount,
            measuredFeatureEdges = measuredFeatureEdgeCount,
            unmeasuredFeatureEdges = unmeasuredFeatureEdgeCount,
            check.names = FALSE
        ),
        assayPairs = S4Vectors::DataFrame(assayPairTable, check.names = FALSE),
        edgeTypes = S4Vectors::DataFrame(edgeTypeTable, check.names = FALSE),
        edgeDirections = S4Vectors::DataFrame(edgeDirectionTable, check.names = FALSE),
        knowledgeSources = S4Vectors::DataFrame(sourceTable, check.names = FALSE),
        evidenceScores = S4Vectors::DataFrame(evidenceSummary, check.names = FALSE),
        warnings = warningMessages
    )
}

#' Read a Knowledge Network from Disk
#'
#' Read a delimited prior-knowledge table, optionally remap the input
#' column names, and standardize it for downstream CorNetto functions.
#'
#' @param filePath Path to a delimited text file.
#' @param columnMapping Optional named character vector mapping standard
#'   CorNetto column names to columns present in the file.
#' @param delimiter Optional delimiter. When `NULL`, the delimiter is
#'   inferred from the file extension.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   standardized network.
#' @param resultName Name used when storing the knowledge network.
#' @param storeResult Whether to store the result in `analysisData`.
#' @param edgeType Default edge type when the column is missing.
#' @param edgeDirection Default edge direction when the column is
#'   missing.
#' @param knowledgeSource Default source label when the column is
#'   missing.
#'
#' @return A standardized `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @examples
#' knowledgePath <- system.file(
#'     "extdata", "example_knowledge_network.csv", package = "CorNetto"
#' )
#' knowledgeNetwork <- readKnowledgeNetwork(knowledgePath)
#' knowledgeNetwork
#' @export
readKnowledgeNetwork <- function(
    filePath,
    columnMapping = NULL,
    delimiter = NULL,
    analysisData = NULL,
    resultName = "knowledgeNetwork",
    storeResult = FALSE,
    edgeType = "association",
    edgeDirection = "undirected",
    knowledgeSource = "userSupplied"
) {
    .assertScalarCharacter(filePath, "filePath")
    if (!file.exists(filePath)) {
        stop("`filePath` does not exist: ", filePath, call. = FALSE)
    }

    delimiter <- .matchDelimiter(filePath, delimiter)
    ## Blank cells are read as NA so that a supplied-but-empty column falls back
    ## to the same default as an absent one in .standardizeKnowledgeNetwork().
    knowledgeNetwork <- utils::read.delim(
        file = filePath,
        sep = delimiter,
        na.strings = c("NA", ""),
        stringsAsFactors = FALSE,
        check.names = FALSE,
        encoding = "UTF-8"
    )
    knowledgeNetwork <- .asPlainDataFrame(knowledgeNetwork)

    if (!is.null(columnMapping)) {
        for (standardName in names(columnMapping)) {
            sourceName <- columnMapping[[standardName]]
            if (sourceName %in% names(knowledgeNetwork)) {
                names(knowledgeNetwork)[
                    names(knowledgeNetwork) == sourceName
                ] <- standardName
            }
        }
    }

    validateKnowledgeNetwork(
        knowledgeNetwork = knowledgeNetwork,
        analysisData = analysisData,
        resultName = resultName,
        storeResult = storeResult,
        edgeType = edgeType,
        edgeDirection = edgeDirection,
        knowledgeSource = knowledgeSource,
        quiet = TRUE
    )
}

#' Combine Multiple Knowledge Networks
#'
#' Row-bind and standardize multiple user-supplied prior-knowledge
#' networks.
#'
#' @param ... Knowledge-network objects or lists of knowledge-network
#'   objects.
#' @param removeDuplicates Whether to remove duplicated interactions.
#' @param analysisData Optional `MultiAssayExperiment` used to store the
#'   combined network.
#' @param resultName Name used when storing the combined network.
#' @param storeResult Whether to store the result in `analysisData`.
#'
#' @return A standardized `DataFrame` or an updated
#'   `MultiAssayExperiment`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' combined <- combineKnowledgeNetworks(
#'     knowledgeNetwork[1:3, ],
#'     knowledgeNetwork[4:6, ]
#' )
#' combined
#' @export
combineKnowledgeNetworks <- function(
    ...,
    removeDuplicates = TRUE,
    analysisData = NULL,
    resultName = "combinedKnowledgeNetwork",
    storeResult = FALSE
) {
    .assertScalarLogical(removeDuplicates, "removeDuplicates")
    .assertScalarLogical(storeResult, "storeResult")
    networkList <- .flattenEdgeTables(...)
    if (!length(networkList)) {
        stop("At least one knowledge network must be supplied.", call. = FALSE)
    }

    standardizedNetworks <- lapply(
        networkList,
        validateKnowledgeNetwork,
        storeResult = FALSE,
        quiet = TRUE
    )
    combinedNetwork <- do.call(rbind, standardizedNetworks)
    combinedNodeTable <- .mergeNodeTables(standardizedNetworks, edgeTable = combinedNetwork)

    if (removeDuplicates) {
        combinedNetwork <- unique(.asPlainDataFrame(combinedNetwork))
        combinedNetwork <- S4Vectors::DataFrame(combinedNetwork, check.names = FALSE)
    }

    combinedNetwork <- .coerceStandardEdgeTable(combinedNetwork, nodeTable = combinedNodeTable)
    if (!storeResult) {
        return(combinedNetwork)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "knowledgeNetworks",
        resultObject = combinedNetwork,
        resultName = resultName
    )
}
