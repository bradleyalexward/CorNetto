.standardizeKnowledgeNetwork <- function(
    knowledgeNetwork,
    edgeType = "association",
    edgeDirection = "undirected",
    knowledgeSource = "userSupplied"
) {
    knowledgeNetwork <- as.data.frame(knowledgeNetwork, stringsAsFactors = FALSE, check.names = FALSE)
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
    if (!"fromFeatureName" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$fromFeatureName <- knowledgeNetwork$fromFeatureIdentifier
    }
    if (!"toFeatureName" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$toFeatureName <- knowledgeNetwork$toFeatureIdentifier
    }
    if (!"isDirected" %in% names(knowledgeNetwork)) {
        knowledgeNetwork$isDirected <-
            !(tolower(knowledgeNetwork$edgeDirection) %in% c("undirected", "association")) &
            !(tolower(knowledgeNetwork$edgeType) %in% c("proteinproteininteraction", "correlation"))
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
        as.numeric(knowledgeNetwork$evidenceScore)
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
#'
#' @return A standardized `DataFrame` or an updated
#'   `MultiAssayExperiment` when `storeResult = TRUE`.
#' @export
validateKnowledgeNetwork <- function(
    knowledgeNetwork,
    analysisData = NULL,
    resultName = "knowledgeNetwork",
    storeResult = FALSE,
    edgeType = "association",
    edgeDirection = "undirected",
    knowledgeSource = "userSupplied"
) {
    standardizedNetwork <- .standardizeKnowledgeNetwork(
        knowledgeNetwork = knowledgeNetwork,
        edgeType = edgeType,
        edgeDirection = edgeDirection,
        knowledgeSource = knowledgeSource
    )

    if (anyNA(standardizedNetwork$fromFeatureIdentifier) || anyNA(standardizedNetwork$toFeatureIdentifier)) {
        stop("Knowledge networks cannot contain missing feature identifiers.", call. = FALSE)
    }
    if (anyNA(standardizedNetwork$fromAssayName) || anyNA(standardizedNetwork$toAssayName)) {
        stop("Knowledge networks cannot contain missing assay names.", call. = FALSE)
    }

    if (!storeResult) {
        return(standardizedNetwork)
    }

    .assertMultiAssayExperiment(analysisData)
    .storeCorNettoResults(
        analysisData = analysisData,
        slotName = "knowledgeNetworks",
        resultObject = standardizedNetwork,
        resultName = resultName
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
    knowledgeNetwork <- readr::read_delim(
      file = filePath,
      delim = delimiter,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal"
    )
    knowledgeNetwork <- as.data.frame(knowledgeNetwork, stringsAsFactors = FALSE)

    if (!is.null(columnMapping)) {
        for (standardName in names(columnMapping)) {
            sourceName <- columnMapping[[standardName]]
            if (sourceName %in% names(knowledgeNetwork)) {
                names(knowledgeNetwork)[names(knowledgeNetwork) == sourceName] <- standardName
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
        knowledgeSource = knowledgeSource
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
#' @export
combineKnowledgeNetworks <- function(
    ...,
    removeDuplicates = TRUE,
    analysisData = NULL,
    resultName = "combinedKnowledgeNetwork",
    storeResult = FALSE
) {
    networkList <- .flattenEdgeTables(...)
    if (!length(networkList)) {
        stop("At least one knowledge network must be supplied.", call. = FALSE)
    }

    standardizedNetworks <- lapply(networkList, validateKnowledgeNetwork, storeResult = FALSE)
    combinedNetwork <- do.call(rbind, standardizedNetworks)
    combinedNodeTable <- .mergeNodeTables(standardizedNetworks, edgeTable = combinedNetwork)

    if (removeDuplicates) {
        combinedNetwork <- unique(as.data.frame(combinedNetwork, stringsAsFactors = FALSE))
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
