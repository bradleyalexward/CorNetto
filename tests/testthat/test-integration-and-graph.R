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

test_that("mirroring undirected edges is idempotent", {
    duplicateUndirected <- getFromNamespace(".duplicateUndirectedEdges", "CorNetto")
    coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")

    network <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("A", "B"),
        toFeatureIdentifier = c("B", "C"),
        fromFeatureName = c("A", "B"),
        toFeatureName = c("B", "C"),
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeType = "correlation",
        isDirected = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    once <- duplicateUndirected(network)
    expect_equal(nrow(once), 4L)

    # Composing combineNetworks(includeReverseEdges = TRUE) with
    # createNetworkGraph() must not keep adding copies.
    twice <- duplicateUndirected(once)
    expect_equal(nrow(twice), 4L)
})

test_that("undirected edges stay symmetric in a mixed directed graph", {
    coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")

    # igraph cannot mix directed and undirected edges, so A-B must be reachable
    # from B even though the graph is directed by C->D.
    network <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("A", "C"),
        toFeatureIdentifier = c("B", "D"),
        fromFeatureName = c("A", "C"),
        toFeatureName = c("B", "D"),
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeType = "correlation",
        isDirected = c(FALSE, TRUE),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    graph <- createNetworkGraph(network)
    expect_true(igraph::is_directed(graph))
    expect_equal(
        length(igraph::neighbors(graph, "protein::B", mode = "out")), 1L
    )
    expect_equal(
        length(igraph::neighbors(graph, "protein::D", mode = "out")), 0L
    )

    unmirrored <- createNetworkGraph(network, mirrorUndirectedEdges = FALSE)
    expect_equal(igraph::gsize(unmirrored), 2L)
    expect_equal(
        length(igraph::neighbors(unmirrored, "protein::B", mode = "out")), 0L
    )
})

test_that("calculateRewiringScores rejects a mirrored network", {
    coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")

    network <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("A", "B"),
        toFeatureIdentifier = c("B", "C"),
        fromFeatureName = c("A", "B"),
        toFeatureName = c("B", "C"),
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeWeight = c(2, 2),
        isDirected = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    expect_silent(calculateRewiringScores(network, storeResult = FALSE))

    mirrored <- combineNetworks(networkList = network, includeReverseEdges = TRUE)
    expect_error(
        calculateRewiringScores(mirrored, storeResult = FALSE),
        "duplicate or mirrored"
    )
})

test_that("base data-frame edge tables are not split into columns", {
    network <- as.data.frame(exampleKnowledgeNetwork())

    combined <- combineNetworks(
        networkList = network,
        includeReverseEdges = FALSE
    )
    knowledge <- combineKnowledgeNetworks(network)

    expect_equal(nrow(combined), nrow(network))
    expect_equal(nrow(knowledge), nrow(network))
    expect_false(anyNA(combined$fromFeatureIdentifier))
})

test_that("node keys reject ambiguous separator-containing components", {
    nodeKey <- getFromNamespace(".nodeKey", "CorNetto")

    expect_error(nodeKey("rna::gene", "TP53"), "cannot contain")
    expect_error(nodeKey("rna", "gene::TP53"), "cannot contain")
})

test_that("graph conversion retains supplied node annotations", {
    network <- exampleKnowledgeNetwork()[1, , drop = FALSE]
    nodeTable <- S4Vectors::DataFrame(
        nodeKey = c("protein::P1", "protein::P2"),
        nodeIdentifier = c("P1", "P2"),
        nodeName = c("P1", "P2"),
        assayName = c("protein", "protein"),
        pathway = c("A", "B"),
        check.names = FALSE
    )

    graph <- createNetworkGraph(network, nodeTable = nodeTable)
    expect_identical(igraph::vertex_attr(graph, "pathway"), c("A", "B"))
})

test_that("public scalar arguments fail clearly", {
    analysisData <- exampleAnalysisData()
    expect_error(
        getCorNettoResult(analysisData, resultType = c("knowledgeNetworks", "correlationResults")),
        "single non-empty character"
    )
    expect_error(
        createCorrelationNetwork(
            analysisData,
            assayName = "protein",
            groupColumn = "clinicalGroup",
            groupLevel = "Recovered",
            sampleIds = "S01",
            storeResult = FALSE
        ),
        "either `sampleIds`"
    )

    rewiringTable <- S4Vectors::DataFrame(
        nodeKey = "protein::P1",
        rootMeanSquareRewiringScore = 1
    )
    expect_error(
        plotRewiringScores(rewiringTable, scoreColumn = c("x", "y")),
        "single non-empty character"
    )
})

test_that("rewiring plots put small permutation tails first", {
    rewiringTable <- S4Vectors::DataFrame(
        nodeKey = paste0("protein::P", 1:3),
        permutationTailProbability = c(0.5, 0.01, 0.05),
        check.names = FALSE
    )
    plotFile <- tempfile(fileext = ".pdf")
    grDevices::pdf(plotFile)
    on.exit(grDevices::dev.off(), add = TRUE)

    plotted <- plotRewiringScores(
        rewiringTable,
        scoreColumn = "permutationTailProbability",
        topN = 2
    )

    expect_identical(as.character(plotted$nodeKey), c("protein::P2", "protein::P3"))
})

test_that("written network tables round-trip in both formats", {
    networkTables <- prepareCytoscapeTables(exampleKnowledgeNetwork())

    for (fileFormat in c("tsv", "csv")) {
        outputDir <- file.path(tempdir(), paste0("cornetto-roundtrip-", fileFormat))
        writtenFiles <- writeNetworkTables(
            networkTables,
            directoryPath = outputDir,
            fileFormat = fileFormat
        )
        expect_true(all(file.exists(writtenFiles)))

        reread <- utils::read.delim(
            writtenFiles[["edges"]],
            sep = if (identical(fileFormat, "csv")) "," else "\t",
            na.strings = c("NA", ""),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        expect_identical(names(reread), names(as.data.frame(networkTables$edges)))
        expect_identical(nrow(reread), nrow(networkTables$edges))
        expect_identical(
            reread$fromFeatureIdentifier,
            as.character(networkTables$edges$fromFeatureIdentifier)
        )
        ## Empty cells, not the literal "NA", so Cytoscape reads them as missing.
        expect_true(all(is.na(reread$correlationValue)))
    }
})
