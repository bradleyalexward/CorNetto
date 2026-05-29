.exampleFilePath <- function(fileName) {
    system.file("extdata", fileName, package = "CorNetto", mustWork = TRUE)
}

.readExampleAssay <- function(fileName) {
    assayTable <- utils::read.csv(
        .exampleFilePath(fileName),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    featureIdentifiers <- as.character(assayTable$featureIdentifier)
    featureNames <- as.character(assayTable$featureName)
    assayMatrix <- as.matrix(
        assayTable[
            ,
            setdiff(names(assayTable), c("featureIdentifier", "featureName")),
            drop = FALSE
        ]
    )
    rownames(assayMatrix) <- featureIdentifiers

    SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = assayMatrix),
        rowData = S4Vectors::DataFrame(
            featureIdentifier = featureIdentifiers,
            featureName = featureNames,
            check.names = FALSE
        )
    )
}

#' Load the Synthetic Example Analysis Data
#'
#' @return A `MultiAssayExperiment`.
#' @examples
#' analysisData <- exampleAnalysisData()
#' names(MultiAssayExperiment::experiments(analysisData))
#' @export
exampleAnalysisData <- function() {
    sampleData <- utils::read.csv(
        .exampleFilePath("example_sample_data.csv"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rownames(sampleData) <- sampleData$sampleId

    assayList <- list(
        protein = .readExampleAssay("example_protein_abundance.csv"),
        transcript = .readExampleAssay("example_transcript_abundance.csv"),
        metabolite = .readExampleAssay("example_metabolite_abundance.csv")
    )

    createAnalysisData(
        assayList = assayList,
        sampleData = sampleData
    )
}

#' Load the Synthetic Example Knowledge Network
#'
#' @return A standardized knowledge-network `DataFrame`.
#' @examples
#' knowledgeNetwork <- exampleKnowledgeNetwork()
#' nrow(knowledgeNetwork)
#' @export
exampleKnowledgeNetwork <- function() {
    readKnowledgeNetwork(
        filePath = .exampleFilePath("example_knowledge_network.csv")
    )
}

#' Load the Synthetic Example Seed Nodes
#'
#' @return A character vector of seed-node identifiers.
#' @examples
#' seedNodes <- exampleSeedNodes()
#' seedNodes
#' @export
exampleSeedNodes <- function() {
    readLines(.exampleFilePath("example_seed_nodes.txt"), warn = FALSE)
}
