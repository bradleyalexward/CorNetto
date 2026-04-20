test_that("example analysis data is a valid MultiAssayExperiment", {
    analysisData <- exampleAnalysisData()
    expect_true(methods::is(analysisData, "MultiAssayExperiment"))

    validated <- validateAnalysisData(analysisData)
    expect_true(methods::is(validated, "MultiAssayExperiment"))

    featureMap <- mapFeatureIdentifiers(validated)
    expect_true(all(c("featureIdentifier", "featureName", "assayName") %in% names(featureMap)))
    expect_gt(nrow(featureMap), 0)
})

test_that("feature filtering reduces assay size", {
    analysisData <- exampleAnalysisData()
    filtered <- filterFeatures(
        analysisData = analysisData,
        assayNames = "protein",
        topVariableFeatures = 3
    )

    proteinMatrix <- MultiAssayExperiment::experiments(filtered)[["protein"]]
    expect_equal(nrow(SummarizedExperiment::assay(proteinMatrix)), 3)
})

test_that("example knowledge network validates and can be combined", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    knowledgeNetwork <- exampleKnowledgeNetwork()
    expect_true(methods::is(knowledgeNetwork, "DataFrame"))
    expect_true(all(standardColumns %in% names(knowledgeNetwork)))
    expect_gt(nrow(getNodeTable(knowledgeNetwork, fallbackToEdges = FALSE)), 0)

    combined <- combineKnowledgeNetworks(knowledgeNetwork, knowledgeNetwork, removeDuplicates = TRUE)
    expect_true(methods::is(combined, "DataFrame"))
    expect_equal(nrow(combined), nrow(knowledgeNetwork))
    expect_gt(nrow(getNodeTable(combined, fallbackToEdges = FALSE)), 0)
})
