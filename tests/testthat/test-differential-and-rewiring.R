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

test_that("differential adjusted p-values are calculated before filtering", {
    analysisData <- exampleAnalysisData()
    resolveGroupSamples <- getFromNamespace(".resolveGroupSamples", "CorNetto")
    resolveFeatureSubset <- getFromNamespace(".resolveFeatureSubset", "CorNetto")
    computeRawDGCA <- getFromNamespace(".computeDifferentialCorrelationDGCA", "CorNetto")
    finalizeResults <- getFromNamespace(".finalizeDifferentialCorrelationResults", "CorNetto")

    groupLevels <- c("PASC", "Recovered")
    group1SampleIds <- resolveGroupSamples(analysisData, "clinicalGroup", groupLevels[[1L]])
    group2SampleIds <- resolveGroupSamples(analysisData, "clinicalGroup", groupLevels[[2L]])
    sharedFeatures <- resolveFeatureSubset(analysisData, "protein")

    rawResults <- computeRawDGCA(
        analysisData = analysisData,
        assayName = "protein",
        sharedFeatures = sharedFeatures,
        group1SampleIds = group1SampleIds,
        group2SampleIds = group2SampleIds,
        groupLevels = groupLevels,
        correlationMethod = "pearson",
        featureNameColumn = "featureName"
    )
    rawData <- as.data.frame(rawResults, stringsAsFactors = FALSE, check.names = FALSE)
    expectedAdjusted <- stats::p.adjust(rawData$pValue, method = "fdr")

    minimumCorrelation <- stats::median(
        abs(c(rawData$group1CorrelationValue, rawData$group2CorrelationValue)),
        na.rm = TRUE
    )
    filteredResults <- finalizeResults(
        edgeTable = rawResults,
        groupLevels = groupLevels,
        minimumAbsoluteCorrelation = minimumCorrelation,
        adjustedPValueThreshold = 1,
        pAdjustMethod = "fdr",
        correlationMethod = "pearson"
    )
    filteredData <- as.data.frame(filteredResults, stringsAsFactors = FALSE, check.names = FALSE)
    expect_gt(nrow(filteredData), 0)

    rawKey <- paste(rawData$fromFeatureIdentifier, rawData$toFeatureIdentifier, sep = "::")
    filteredKey <- paste(filteredData$fromFeatureIdentifier, filteredData$toFeatureIdentifier, sep = "::")
    expect_equal(
        filteredData$adjustedPValue,
        expectedAdjusted[match(filteredKey, rawKey)]
    )
})

test_that("rewiring permutation validation returns empirical node p-values", {
    analysisData <- exampleAnalysisData()

    validationResult <- permuteDifferentialRewiring(
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
    )

    expect_true(all(c(
        "differentialCorrelationTable",
        "differentialCorrelationNetwork",
        "rewiringTable",
        "permutationScores"
    ) %in% names(validationResult)))
    expect_true(all(c(
        "empiricalPValue",
        "empiricalAdjustedPValue",
        "nullMeanRawRewiringScore",
        "nPermutations"
    ) %in% names(validationResult$rewiringTable)))
    expect_true(all(validationResult$rewiringTable$nPermutations == 2))
    expect_true(methods::is(validationResult$permutationScores, "DataFrame"))
})
