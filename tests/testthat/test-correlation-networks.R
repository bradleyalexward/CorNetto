test_that("dense correlation networks support all descriptive methods", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()

    pearsonNetwork <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "pearson",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    spearmanNetwork <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "spearman",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    kendallNetwork <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        correlationMethod = "kendall",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))

    expect_true(all(standardColumns %in% names(pearsonNetwork)))
    expect_true(all(standardColumns %in% names(spearmanNetwork)))
    expect_true(all(standardColumns %in% names(kendallNetwork)))
    expect_true(all(pearsonNetwork$correlationMethod == "pearson"))
    expect_true(all(spearmanNetwork$correlationMethod == "spearman"))
    expect_true(all(kendallNetwork$correlationMethod == "kendall"))
    expect_true(all(c("sampleCount", "group1SampleCount", "group2SampleCount") %in% names(pearsonNetwork)))
    expect_true(all(!is.na(pearsonNetwork$sampleCount)))
    expect_gt(nrow(getNodeTable(pearsonNetwork, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(spearmanNetwork, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(kendallNetwork, fallbackToEdges = FALSE)), 0)
})

test_that("sparse candidate selection supports all, within, and cross scopes", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    selectCandidateEdges <- getFromNamespace(".selectSparseCandidateEdges", "CorNetto")
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))
    measuredKnowledgeNetwork <- knowledgeNetwork[
        knowledgeNetwork$fromAssayName %in% measuredAssays &
            knowledgeNetwork$toAssayName %in% measuredAssays,
        ,
        drop = FALSE
    ]

    sparseAll <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "all"
    )
    sparseWithin <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "withinOmic"
    )
    sparseCross <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        correlationScope = "crossOmic"
    )

    expect_true(all(standardColumns %in% names(sparseAll)))
    expect_true(all(standardColumns %in% names(sparseWithin)))
    expect_true(all(standardColumns %in% names(sparseCross)))
    expect_true(all(sparseWithin$fromAssayName == sparseWithin$toAssayName))
    expect_true(all(sparseCross$fromAssayName != sparseCross$toAssayName))
    expect_true(any(sparseAll$fromAssayName == sparseAll$toAssayName))
    expect_true(any(sparseAll$fromAssayName != sparseAll$toAssayName))
    sparseProteinOnly <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        assayNames = "protein",
        correlationScope = "all"
    )
    expect_true(all(sparseProteinOnly$fromAssayName == "protein"))
    expect_true(all(sparseProteinOnly$toAssayName == "protein"))

    sparseProteinEither <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = measuredKnowledgeNetwork,
        assayNames = "protein",
        assayFilterMode = "either",
        correlationScope = "all"
    )
    expect_true(any(sparseProteinEither$fromAssayName != "protein" | sparseProteinEither$toAssayName != "protein"))
    expect_gt(nrow(getNodeTable(sparseAll, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(sparseWithin, fallbackToEdges = FALSE)), 0)
    expect_gt(nrow(getNodeTable(sparseCross, fallbackToEdges = FALSE)), 0)
})

test_that("combineNetworks can merge group-specific correlation outputs", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()

    recovered <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    pasc <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "PASC",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))

    combined <- combineNetworks(
        correlationNetwork = list(recovered, pasc),
        includeReverseEdges = FALSE
    )
    expect_true(methods::is(combined, "DataFrame"))
    expect_true(all(standardColumns %in% names(combined)))
    expect_gte(nrow(combined), max(nrow(recovered), nrow(pasc)))
    expect_gt(nrow(getNodeTable(combined, fallbackToEdges = FALSE)), 0)
})

test_that("sparse candidate deduplication retains direction and edge type", {
    analysisData <- exampleAnalysisData()
    selectCandidateEdges <- getFromNamespace(".selectSparseCandidateEdges", "CorNetto")
    prior <- S4Vectors::DataFrame(
        fromFeatureIdentifier = c("P1", "P1", "P1"),
        toFeatureIdentifier = c("P2", "P2", "P2"),
        fromAssayName = c("protein", "protein", "protein"),
        toAssayName = c("protein", "protein", "protein"),
        edgeType = c("activation", "inhibition", "activation"),
        edgeDirection = c("directed", "directed", "undirected"),
        isDirected = c(TRUE, TRUE, FALSE),
        check.names = FALSE
    )

    sparseEdges <- selectCandidateEdges(
        analysisData = analysisData,
        knowledgeNetwork = prior,
        correlationScope = "withinOmic"
    )

    expect_equal(nrow(sparseEdges), 3)
})
