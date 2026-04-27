test_that("dense correlation networks support all descriptive methods", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()

    pearsonNetwork <- createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "pearson",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )
    spearmanNetwork <- createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "spearman",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )
    kendallNetwork <- createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "kendall",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )

    expect_true(all(standardColumns %in% names(pearsonNetwork)))
    expect_true(all(standardColumns %in% names(spearmanNetwork)))
    expect_true(all(standardColumns %in% names(kendallNetwork)))
    expect_true(all(pearsonNetwork$correlationMethod == "pearson"))
    expect_true(all(spearmanNetwork$correlationMethod == "spearman"))
    expect_true(all(kendallNetwork$correlationMethod == "kendall"))
    expect_gt(nrow(getNodeTable(pearsonNetwork, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(spearmanNetwork, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(kendallNetwork, fallbackToEdges = FALSE)), 0)
})

test_that("sparse multi-omic correlations support all, within, and cross scopes", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))
    measuredKnowledgeNetwork <- knowledgeNetwork[
        knowledgeNetwork$fromAssayName %in% measuredAssays &
            knowledgeNetwork$toAssayName %in% measuredAssays,
        ,
        drop = FALSE
    ]

    sparseAll <- createSparseMultiOmicCorrelations(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "all",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )
    sparseWithin <- createSparseMultiOmicCorrelations(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "withinOmic",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )
    sparseCross <- createSparseMultiOmicCorrelations(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "crossOmic",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )

    expect_true(all(standardColumns %in% names(sparseAll)))
    expect_true(all(standardColumns %in% names(sparseWithin)))
    expect_true(all(standardColumns %in% names(sparseCross)))
    expect_true(any(sparseAll$correlationScope == "sparseWithinOmic"))
    expect_true(any(sparseAll$correlationScope == "sparseCrossOmic"))
    expect_true(all(sparseWithin$fromAssayName == sparseWithin$toAssayName))
    expect_true(all(sparseCross$fromAssayName != sparseCross$toAssayName))
    expect_gt(nrow(getNodeTable(sparseAll, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(sparseWithin, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(sparseCross, fallbackToEdges = FALSE)), 0)
})

test_that("combined correlation networks can merge dense and sparse outputs", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))
    measuredKnowledgeNetwork <- knowledgeNetwork[
        knowledgeNetwork$fromAssayName %in% measuredAssays &
            knowledgeNetwork$toAssayName %in% measuredAssays,
        ,
        drop = FALSE
    ]

    dense <- createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )
    sparseMulti <- createSparseMultiOmicCorrelations(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "all",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    )

    combined <- combineCorrelationNetworks(dense, sparseMulti)
    expect_true(methods::is(combined, "DataFrame"))
    expect_true(all(standardColumns %in% names(combined)))
    expect_gte(nrow(combined), max(nrow(dense), nrow(sparseMulti)))
    expect_gt(nrow(getNodeTable(combined, fallbackToEdges = FALSE)), 0)
})

test_that("sparse edge permutation validation returns empirical p-values", {
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()

    validationTable <- validateSparseCorrelationEdges(
        analysisData = analysisData,
        knowledgeNetwork = knowledgeNetwork,
        correlationScope = "all",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        nPermutations = 3,
        seed = 1
    )

    expect_true(methods::is(validationTable, "DataFrame"))
    expect_true(all(c(
        "observedCorrelationValue",
        "observedAdjustedPValue",
        "empiricalPValue",
        "empiricalAdjustedPValue",
        "nPermutations",
        "nullModel"
    ) %in% names(validationTable)))
    expect_true(all(validationTable$nPermutations == 3))
})
