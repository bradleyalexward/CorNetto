test_that("integrated network and neighborhood workflow are available", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    setNodeTable <- getFromNamespace(".setStoredNodeTable", "CorNetto")
    coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()

    dense <- suppressWarnings(createCorrelationNetwork(
        analysisData = analysisData,
        assayName = "protein",
        groupColumn = "clinicalGroup",
        groupLevel = "Recovered",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    differentialResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("Recovered", "PASC"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    ))
    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialResults,
        differenceAdjustedPValueThreshold = 1,
        minimumAbsoluteCorrelation = 0
    )

    integrated <- combineNetworks(
        knowledgeNetwork = knowledgeNetwork,
        correlationNetwork = dense,
        differentialCorrelationNetwork = differentialNetwork,
        includeReverseEdges = TRUE
    )
    expect_true(all(standardColumns %in% names(integrated)))
    expect_gt(nrow(integrated), 0)
    expect_gt(nrow(getNodeTable(integrated, fallbackToEdges = FALSE)), 0)

    pathwayNetwork <- createFocusedNetwork(
        networkEdgeTable = integrated,
        seedNodes = exampleSeedNodes(),
        neighborhoodOrder = 1
    )
    expect_true(methods::is(pathwayNetwork, "DataFrame"))
    expect_gt(nrow(getNodeTable(pathwayNetwork, fallbackToEdges = FALSE)), 0)

    cytoscapeTables <- prepareCytoscapeTables(pathwayNetwork)
    expect_true(all(c("nodes", "edges") %in% names(cytoscapeTables)))

    graphObject <- createNetworkGraph(pathwayNetwork)
    expect_true(inherits(graphObject, "igraph"))

    isolateNetwork <- coerceEdgeTable(
        data.frame(
            fromFeatureIdentifier = "seedA",
            toFeatureIdentifier = "seedB",
            fromFeatureName = "seedA",
            toFeatureName = "seedB",
            fromAssayName = "RNA",
            toAssayName = "RNA",
            edgeType = "correlation",
            edgeDirection = "positive",
            sourceType = "correlation",
            correlationScope = "withinOmic",
            correlationMethod = "pearson",
            edgeWeight = 1,
            isDirected = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    )
    isolateNetwork <- setNodeTable(
        isolateNetwork,
        S4Vectors::DataFrame(
            nodeKey = c("RNA::seedA", "RNA::seedB", "RNA::seedC"),
            nodeIdentifier = c("seedA", "seedB", "seedC"),
            nodeName = c("seedA", "seedB", "seedC"),
            assayName = c("RNA", "RNA", "RNA"),
            check.names = FALSE
        )
    )

    focusedWithIsolate <- createFocusedNetwork(
        networkEdgeTable = isolateNetwork,
        seedNodes = "seedC",
        neighborhoodOrder = 0,
        dropIsolatedNodes = FALSE
    )
    focusedNodeTable <- getNodeTable(focusedWithIsolate, fallbackToEdges = FALSE)
    expect_equal(nrow(focusedWithIsolate), 0)
    expect_equal(nrow(focusedNodeTable), 1)
    expect_equal(focusedNodeTable$nodeIdentifier, "seedC")

    focusedCytoscape <- prepareCytoscapeTables(focusedWithIsolate)
    expect_equal(nrow(focusedCytoscape$nodes), 1)
    expect_equal(nrow(focusedCytoscape$edges), 0)
    rewiringTable <- S4Vectors::DataFrame(
        nodeKey = "RNA::seedC",
        nodeIdentifier = "seedC",
        assayName = "RNA",
        totalConnections = 0L,
        rawRewiringScore = 0,
        check.names = FALSE
    )
    focusedCytoscapeWithRewiring <- prepareCytoscapeTables(
        networkEdgeTable = focusedWithIsolate,
        rewiringTable = rewiringTable
    )
    expect_true(all(c("nodeIdentifier", "nodeName", "assayName") %in% names(focusedCytoscapeWithRewiring$nodes)))
    expect_false(any(c("nodeIdentifier.x", "assayName.x") %in% names(focusedCytoscapeWithRewiring$nodes)))
    expect_true("rawRewiringScore" %in% names(focusedCytoscapeWithRewiring$nodes))

    isolateGraph <- createNetworkGraph(focusedWithIsolate)
    expect_true(inherits(isolateGraph, "igraph"))
    expect_equal(igraph::gorder(isolateGraph), 1)
    expect_equal(igraph::gsize(isolateGraph), 0)

    focusedWithoutIsolate <- createFocusedNetwork(
        networkEdgeTable = isolateNetwork,
        seedNodes = "seedC",
        neighborhoodOrder = 0,
        dropIsolatedNodes = TRUE
    )
    expect_equal(nrow(getNodeTable(focusedWithoutIsolate, fallbackToEdges = FALSE)), 0)

    focusedByKey <- createFocusedNetwork(
        networkEdgeTable = isolateNetwork,
        seedNodes = "RNA::seedC",
        nodeIdType = "key",
        neighborhoodOrder = 0,
        dropIsolatedNodes = FALSE
    )
    expect_equal(nrow(getNodeTable(focusedByKey, fallbackToEdges = FALSE)), 1)

    filteredByKey <- filterNetworkByNodes(
        networkEdgeTable = isolateNetwork,
        nodes = "RNA::seedA",
        nodeIdType = "key",
        mode = "either"
    )
    expect_equal(nrow(filteredByKey), 1)
})
