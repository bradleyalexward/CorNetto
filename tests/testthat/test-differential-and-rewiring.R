test_that("differential correlation and rewiring workflow returns expected columns", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()

    differentialResults <- testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    )

    expect_true(all(standardColumns %in% names(differentialResults)))
    expect_gt(nrow(differentialResults), 0)
    expect_gt(nrow(getNodeTable(differentialResults, fallbackToEdges = FALSE)), 0)
    expect_true(all(!is.na(differentialResults$fromFeatureIdentifier)))
    expect_true(all(!is.na(differentialResults$toFeatureIdentifier)))

    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = 1,
        minimumAbsoluteCorrelation = 0
    )
    expect_true(all(standardColumns %in% names(differentialNetwork)))
    expect_true(all(!is.na(differentialNetwork$edgeDirection)))
    expect_gt(nrow(getNodeTable(differentialNetwork, fallbackToEdges = FALSE)), 0)

    rewiringTable <- calculateRewiringScores(
        differentialCorrelationNetwork = differentialNetwork,
        storeResult = FALSE
    )
    resolveRewiringPlotLabels <- getFromNamespace(".resolveRewiringPlotLabels", "CorNetto")
    expect_true(all(c(
        "nodeKey",
        "nodeIdentifier",
        "nodeName",
        "assayName",
        "totalConnections",
        "rawRewiringScore",
        "rootMeanSquareRewiringScore",
        "degreeMatchedZScore",
        "degreeMatchedScaledScore"
    ) %in% names(rewiringTable)))
    expect_gt(nrow(rewiringTable), 0)
    expect_true(all(!is.na(rewiringTable$nodeName)))
    expect_identical(
        resolveRewiringPlotLabels(rewiringTable, label = "name"),
        as.character(rewiringTable$nodeName)
    )
    expect_identical(
        resolveRewiringPlotLabels(rewiringTable, label = "identifier"),
        as.character(rewiringTable$nodeIdentifier)
    )
})

test_that("differential correlation handles empty subsets and supports spearman", {
    analysisData <- exampleAnalysisData()

    emptyResults <- testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "metabolite",
        featureSubset = character(),
        storeResult = FALSE
    )
    expect_equal(nrow(emptyResults), 0)

    spearmanResults <- testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        correlationMethod = "spearman",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    )

    expect_gt(nrow(spearmanResults), 0)
    expect_true(all(spearmanResults$correlationMethod == "spearman"))
    expect_true(all(!is.na(spearmanResults$fromFeatureIdentifier)))
    expect_true(all(!is.na(spearmanResults$toFeatureIdentifier)))
})
