.exampleFilePath <- function(fileName) {
    system.file("extdata", fileName, package = "CorNetto", mustWork = TRUE)
}

#' Load the Synthetic Example Analysis Data
#'
#' @return A `MultiAssayExperiment`.
#' @export
exampleAnalysisData <- function() {
    sampleData <- utils::read.csv(
        .exampleFilePath("example_sample_data.csv"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rownames(sampleData) <- sampleData$sampleId

    assayList <- list(
        protein = readAssayData(
            filePath = .exampleFilePath("example_protein_abundance.csv"),
            featureIdentifierColumn = "featureIdentifier",
            featureNameColumn = "featureName"
        ),
        transcript = readAssayData(
            filePath = .exampleFilePath("example_transcript_abundance.csv"),
            featureIdentifierColumn = "featureIdentifier",
            featureNameColumn = "featureName"
        ),
        metabolite = readAssayData(
            filePath = .exampleFilePath("example_metabolite_abundance.csv"),
            featureIdentifierColumn = "featureIdentifier",
            featureNameColumn = "featureName"
        )
    )

    createAnalysisData(
        assayList = assayList,
        sampleData = sampleData
    )
}

#' Load the Synthetic Example Knowledge Network
#'
#' @return A standardized knowledge-network `DataFrame`.
#' @export
exampleKnowledgeNetwork <- function() {
    readKnowledgeNetwork(
        filePath = .exampleFilePath("example_knowledge_network.csv")
    )
}

#' Load the Synthetic Example Seed Nodes
#'
#' @return A character vector of seed-node identifiers.
#' @export
exampleSeedNodes <- function() {
    readLines(.exampleFilePath("example_seed_nodes.txt"), warn = FALSE)
}
