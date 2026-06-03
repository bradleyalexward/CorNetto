test_that("differential correlation and rewiring workflow returns expected columns", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    analysisData <- exampleAnalysisData()

    differentialResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    ))

    expect_true(all(standardColumns %in% names(differentialResults)))
    expect_gt(nrow(differentialResults), 0)
    expect_gt(nrow(getNodeTable(differentialResults, fallbackToEdges = FALSE)), 0)
    expect_true(all(!is.na(differentialResults$fromFeatureIdentifier)))
    expect_true(all(!is.na(differentialResults$toFeatureIdentifier)))
    expect_true(all(c("group1SampleCount", "group2SampleCount") %in% names(differentialResults)))

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
        "degreeMatchedZScore"
    ) %in% names(rewiringTable)))
    expect_false("degreeMatchedScaledScore" %in% names(rewiringTable))
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

    emptyResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "metabolite",
        featureSubset = character(),
        storeResult = FALSE
    ))
    expect_equal(nrow(emptyResults), 0)

    spearmanResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        correlationMethod = "spearman",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    ))

    expect_gt(nrow(spearmanResults), 0)
    expect_true(all(spearmanResults$correlationMethod == "spearman"))
    expect_true(all(!is.na(spearmanResults$fromFeatureIdentifier)))
    expect_true(all(!is.na(spearmanResults$toFeatureIdentifier)))
})

test_that("candidate-edge differential correlation ignores unmeasured assays", {
    analysisData <- exampleAnalysisData()
    knowledgeNetwork <- exampleKnowledgeNetwork()
    measuredAssays <- names(MultiAssayExperiment::experiments(analysisData))

    differentialResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = knowledgeNetwork,
        correlationMethod = "pearson",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    ))

    expect_gt(nrow(differentialResults), 0)
    expect_true(all(differentialResults$fromAssayName %in% measuredAssays))
    expect_true(all(differentialResults$toAssayName %in% measuredAssays))
})

test_that("rewiring-first filtering couples strength and significance within group", {
    finalizeResults <- getFromNamespace(".finalizeDifferentialCorrelationResults", "CorNetto")

    rawResults <- S4Vectors::DataFrame(
        fromFeatureIdentifier = c("A", "B", "C", "D"),
        toFeatureIdentifier = c("E", "F", "G", "H"),
        fromFeatureName = c("A", "B", "C", "D"),
        toFeatureName = c("E", "F", "G", "H"),
        fromAssayName = rep("protein", 4),
        toAssayName = rep("protein", 4),
        group1CorrelationValue = c(0.5, 0.5, 0.1, 0.1),
        group2CorrelationValue = c(0.1, 0.1, 0.5, 0.1),
        group1PValue = c(0.001, 0.9, 0.9, 0.001),
        group2PValue = c(0.9, 0.001, 0.001, 0.001),
        pValue = c(0.4, 0.2, 0.1, 0.001),
        zScoreDifference = c(1, 2, 3, 4),
        check.names = FALSE
    )

    filteredResults <- finalizeResults(
        edgeTable = rawResults,
        groupLevels = c("Recovered", "PASC"),
        minimumAbsoluteCorrelation = 0.3,
        adjustedPValueThreshold = 0.01,
        pAdjustMethod = "fdr",
        correlationMethod = "pearson"
    )
    filteredData <- as.data.frame(filteredResults)
    prefilteredData <- as.data.frame(rawResults)[1:3, , drop = FALSE]

    expect_equal(filteredData$fromFeatureIdentifier, c("A", "C"))
    expect_false("B" %in% filteredData$fromFeatureIdentifier)
    expect_false("D" %in% filteredData$fromFeatureIdentifier)
    expect_equal(
        filteredData$adjustedPValue,
        stats::p.adjust(prefilteredData$pValue, method = "fdr")[c(1, 3)]
    )
    expect_equal(
        filteredData$group1AdjustedPValue[[1]],
        stats::p.adjust(prefilteredData$group1PValue, method = "fdr")[[1]]
    )
    expect_equal(
        filteredData$group2AdjustedPValue[[2]],
        stats::p.adjust(prefilteredData$group2PValue, method = "fdr")[[3]]
    )
})

test_that("weighted networks default to rewiring-first signed z-score weights", {
    coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")
    analysisData <- exampleAnalysisData()
    differentialResults <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("A", "B"),
        toFeatureIdentifier = c("C", "D"),
        fromFeatureName = c("A", "B"),
        toFeatureName = c("C", "D"),
        fromAssayName = rep("protein", 2),
        toAssayName = rep("protein", 2),
        group1CorrelationValue = c(0.5, 0.5),
        group2CorrelationValue = c(0.1, -0.5),
        adjustedPValue = c(0.8, 0.9),
        zScoreDifference = c(-2, 3),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        minimumAbsoluteCorrelation = 0.3
    )
    expect_equal(nrow(differentialNetwork), 2)
    expect_equal(differentialNetwork$edgeWeight, c(-2, 3))

    filteredNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = 0.05,
        minimumAbsoluteCorrelation = 0.3
    )
    expect_equal(nrow(filteredNetwork), 0)

    updated <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = 0.05,
        analysisData = analysisData,
        storeResult = TRUE
    )
    expect_true(methods::is(updated, "MultiAssayExperiment"))
    expect_equal(
        nrow(differentialCorrelationNetworks(updated)$differentialCorrelationNetwork),
        0
    )
})

test_that("rewiring permutation validation returns empirical node p-values", {
    analysisData <- exampleAnalysisData()

    validationResult <- suppressWarnings(permuteRewiringScores(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        differenceAdjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        nPermutations = 2,
        seed = 1,
        keepPermutationScores = TRUE
    ))

    expect_true(all(c(
        "differentialCorrelationTable",
        "differentialCorrelationNetwork",
        "rewiringTable",
        "permutationScores"
    ) %in% names(validationResult)))
    expect_true(all(c(
        "empiricalPValue",
        "empiricalAdjustedPValue",
        "nullMeanScore",
        "scoreColumn",
        "nPermutations"
    ) %in% names(validationResult$rewiringTable)))
    expect_true(all(validationResult$rewiringTable$nPermutations == 2))
    expect_true(methods::is(validationResult$permutationScores, "DataFrame"))
    expect_false("degreeMatchedScaledScore" %in% names(validationResult$permutationScores))
})

test_that("weighted differential networks can be stored separately from raw results", {
    analysisData <- exampleAnalysisData()
    differentialResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        storeResult = FALSE
    ))

    updated <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = 1,
        analysisData = analysisData,
        storeResult = TRUE
    )

    storedNetworks <- differentialCorrelationNetworks(updated)
    expect_true("differentialCorrelationNetwork" %in% names(storedNetworks))
    expect_true(all(storedNetworks$differentialCorrelationNetwork$sourceType == "differentialCorrelation"))
    expect_identical(
        getCorNettoResult(
            updated,
            "differentialCorrelationNetworks",
            "differentialCorrelationNetwork"
        ),
        storedNetworks$differentialCorrelationNetwork
    )
})
