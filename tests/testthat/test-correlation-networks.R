test_that("dense correlation networks support all descriptive methods", {
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

test_that("combineNetworks can merge group-specific correlation outputs", {
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
