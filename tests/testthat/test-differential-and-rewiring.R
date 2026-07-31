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

test_that("NULL assayName runs dense within-omic testing across all assays", {
    analysisData <- exampleAnalysisData()
    differentialResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))

    expect_equal(nrow(differentialResults), 26)
    expect_equal(sum(differentialResults$fromAssayName == "protein"), 10)
    expect_equal(sum(differentialResults$fromAssayName == "transcript"), 10)
    expect_equal(sum(differentialResults$fromAssayName == "metabolite"), 6)
    expect_true(all(
        differentialResults$fromAssayName == differentialResults$toAssayName
    ))
    expect_true(all(differentialResults$correlationScope == "withinOmic"))
    expect_equal(
        differentialResults$adjustedPValue,
        stats::p.adjust(differentialResults$pValue, method = "fdr")
    )
    expect_equal(
        differentialResults$group1AdjustedPValue,
        stats::p.adjust(differentialResults$group1PValue, method = "fdr")
    )
    expect_equal(
        differentialResults$group2AdjustedPValue,
        stats::p.adjust(differentialResults$group2PValue, method = "fdr")
    )

    differentialNetwork <- createDifferentialCorrelationNetwork(
        differentialCorrelationTable = differentialResults,
        differenceAdjustedPValueThreshold = 1,
        minimumAbsoluteCorrelation = 0
    )
    rewiringTable <- calculateRewiringScores(differentialNetwork)
    expect_equal(nrow(differentialNetwork), 26)
    expect_equal(nrow(rewiringTable), 14)

    subsetResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        featureSubset = c("P1", "P2", "T1", "T2", "M1", "M2"),
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    expect_equal(nrow(subsetResults), 3)
    expect_true(all(subsetResults$fromAssayName == subsetResults$toAssayName))

    filteredAnalysis <- suppressWarnings(filterFeatures(
        analysisData,
        featureSubset = c("P1", "P2", "T1", "M1", "M2")
    ))
    filteredResults <- suppressWarnings(testDifferentialCorrelation(
        analysisData = filteredAnalysis,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))
    expect_equal(nrow(filteredResults), 2)
    expect_false(any(filteredResults$fromAssayName == "transcript"))
    expect_true(all(filteredResults$fromAssayName == filteredResults$toAssayName))
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

test_that("candidate-edge differential correlation drops self-loops", {
    analysisData <- exampleAnalysisData()
    candidateEdges <- S4Vectors::DataFrame(
        fromFeatureIdentifier = c("P1", "P1"),
        toFeatureIdentifier = c("P1", "P2"),
        fromAssayName = c("protein", "protein"),
        toAssayName = c("protein", "protein"),
        check.names = FALSE
    )

    expect_warning(
        differentialResults <- testDifferentialCorrelation(
            analysisData = analysisData,
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            candidateEdgeTable = candidateEdges,
            minimumAbsoluteCorrelation = 0,
            adjustedPValueThreshold = 1,
            storeResult = FALSE
        ),
        "Dropping 1 self-loop candidate edge"
    )

    expect_equal(nrow(differentialResults), 1)
    expect_false(any(
        differentialResults$fromFeatureIdentifier ==
            differentialResults$toFeatureIdentifier &
            differentialResults$fromAssayName == differentialResults$toAssayName
    ))
})

test_that("p-values are adjusted before correlation-strength filtering", {
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
    rawData <- as.data.frame(rawResults)

    expect_equal(filteredData$fromFeatureIdentifier, c("A", "C"))
    expect_false("B" %in% filteredData$fromFeatureIdentifier)
    expect_false("D" %in% filteredData$fromFeatureIdentifier)
    expect_equal(
        filteredData$adjustedPValue,
        stats::p.adjust(rawData$pValue, method = "fdr")[c(1, 3)]
    )
    expect_equal(
        filteredData$group1AdjustedPValue[[1]],
        stats::p.adjust(rawData$group1PValue, method = "fdr")[[1]]
    )
    expect_equal(
        filteredData$group2AdjustedPValue[[2]],
        stats::p.adjust(rawData$group2PValue, method = "fdr")[[3]]
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

test_that("all-pairs permutation validation returns conditional rankings", {
    analysisData <- exampleAnalysisData()

    expect_warning(
        validationResult <- permuteRewiringScores(
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
        ),
        "prespecified `candidateEdgeTable`"
    )

    expect_true(all(c(
        "differentialCorrelationTable",
        "differentialCorrelationNetwork",
        "rewiringTable",
        "inferenceStatus",
        "permutationScores"
    ) %in% names(validationResult)))
    expect_true(all(c(
        "permutationTailProbability",
        "adjustedPermutationTailProbability",
        "nullMeanScore",
        "scoreColumn",
        "nPermutations",
        "inferenceStatus"
    ) %in% names(validationResult$rewiringTable)))
    expect_identical(validationResult$inferenceStatus, "conditional ranking")
    expect_true(all(validationResult$rewiringTable$inferenceStatus == "conditional ranking"))
    expect_false(any(c("empiricalPValue", "empiricalAdjustedPValue") %in%
        names(validationResult$rewiringTable)))
    expect_true(all(validationResult$rewiringTable$nPermutations == 2))
    expect_true(methods::is(validationResult$permutationScores, "DataFrame"))
    expect_false("degreeMatchedScaledScore" %in% names(validationResult$permutationScores))
})

test_that("all-assay dense permutation validation stays within omics", {
    expect_warning(
        validationResult <- permuteRewiringScores(
            analysisData = exampleAnalysisData(),
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            minimumAbsoluteCorrelation = 0,
            adjustedPValueThreshold = 1,
            differenceAdjustedPValueThreshold = 1,
            nPermutations = 2,
            seed = 1,
            keepPermutationScores = TRUE
        ),
        "prespecified `candidateEdgeTable`"
    )

    expect_equal(nrow(validationResult$differentialCorrelationTable), 26)
    expect_equal(nrow(validationResult$differentialCorrelationNetwork), 26)
    expect_equal(nrow(validationResult$rewiringTable), 14)
    expect_equal(nrow(validationResult$permutationScores), 28)
    expect_true(all(
        validationResult$differentialCorrelationTable$fromAssayName ==
            validationResult$differentialCorrelationTable$toAssayName
    ))
    expect_true(all(
        validationResult$differentialCorrelationNetwork$fromAssayName ==
            validationResult$differentialCorrelationNetwork$toAssayName
    ))
    expect_identical(validationResult$inferenceStatus, "conditional ranking")
})

test_that("a complete fixed candidate universe supports randomization inference", {
    analysisData <- exampleAnalysisData()
    candidateEdge <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )

    validationResult <- suppressWarnings(permuteRewiringScores(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = candidateEdge,
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        differenceAdjustedPValueThreshold = NULL,
        pAdjustMethod = "fdr",
        nPermutations = 3,
        seed = 1
    ))

    rewiringTable <- as.data.frame(validationResult$rewiringTable)
    expect_identical(validationResult$inferenceStatus, "randomization p-value")
    expect_true(all(rewiringTable$inferenceStatus == "randomization p-value"))
    expect_true(all(rewiringTable$permutationTailProbability >= 1 / 4))
    expect_true(all(rewiringTable$permutationTailProbability <= 1))
    expect_equal(
        rewiringTable$adjustedPermutationTailProbability,
        stats::p.adjust(rewiringTable$permutationTailProbability, method = "fdr")
    )

    storedData <- suppressWarnings(permuteRewiringScores(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = candidateEdge,
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        nPermutations = 1,
        seed = 1,
        resultName = "fixedCandidate",
        storeResult = TRUE
    ))
    storedResult <- validationResults(storedData)[["fixedCandidate"]]
    expect_identical(storedResult$inferenceStatus, "randomization p-value")
})

test_that("unestimable permutation scores are not labelled as p-values", {
    validationResult <- suppressWarnings(permuteRewiringScores(
        analysisData = exampleAnalysisData(),
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = exampleKnowledgeNetwork(),
        correlationMethod = "pearson",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        differenceAdjustedPValueThreshold = NULL,
        scoreColumn = "degreeMatchedZScore",
        nPermutations = 2,
        seed = 1
    ))

    rewiringTable <- as.data.frame(validationResult$rewiringTable)
    incomplete <- is.na(rewiringTable$degreeMatchedZScore) |
        rewiringTable$contributingPermutations < rewiringTable$nPermutations
    expect_true(any(incomplete))
    expect_true(all(
        rewiringTable$inferenceStatus[incomplete] == "conditional ranking"
    ))
    expect_identical(validationResult$inferenceStatus, "conditional ranking")
})

test_that("candidate selection and empty universes have explicit status", {
    analysisData <- exampleAnalysisData()
    candidateEdge <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )

    expect_warning(
        selectedResult <- permuteRewiringScores(
            analysisData = analysisData,
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            candidateEdgeTable = candidateEdge,
            minimumAbsoluteCorrelation = 0.3,
            adjustedPValueThreshold = NULL,
            nPermutations = 1,
            seed = 1
        ),
        "conditional rankings"
    )
    expect_identical(selectedResult$inferenceStatus, "conditional ranking")

    selfLoop <- candidateEdge
    selfLoop$toFeatureIdentifier <- selfLoop$fromFeatureIdentifier
    expect_warning(
        emptyResult <- permuteRewiringScores(
            analysisData = analysisData,
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            candidateEdgeTable = selfLoop,
            minimumAbsoluteCorrelation = 0,
            adjustedPValueThreshold = NULL,
            nPermutations = 1,
            seed = 1
        ),
        "Dropping 1 self-loop"
    )
    expect_identical(emptyResult$inferenceStatus, "conditional ranking")
    expect_equal(nrow(emptyResult$rewiringTable), 0)
    expect_true(all(c(
        "permutationTailProbability",
        "adjustedPermutationTailProbability",
        "inferenceStatus"
    ) %in% names(emptyResult$rewiringTable)))
})

test_that("permutation jobs retain every index in order", {
    makeJobs <- getFromNamespace(".makePermutationJobs", "CorNetto")
    labels <- as.list(seq_len(7L))
    jobs <- makeJobs(
        indices = seq_len(7L),
        permutedLabels = labels,
        nWorkers = 3L
    )

    expect_identical(
        unlist(lapply(jobs, `[[`, "indices"), use.names = FALSE),
        seq_len(7L)
    )
    expect_identical(
        vapply(jobs, function(job) length(job$indices), integer(1L)),
        c(2L, 2L, 3L)
    )
})

test_that("parallel permutation scoring matches serial scoring", {
    analysisData <- exampleAnalysisData()
    candidateEdge <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )
    arguments <- list(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = candidateEdge,
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        nPermutations = 3L,
        seed = 17L,
        keepPermutationScores = TRUE,
        progressEvery = 2L
    )

    set.seed(101)
    rngBefore <- .Random.seed
    serialResult <- do.call(
        permuteRewiringScores,
        c(arguments, list(BPPARAM = BiocParallel::SerialParam()))
    )
    expect_identical(.Random.seed, rngBefore)

    snow <- BiocParallel::SnowParam(
        workers = 2L,
        type = "SOCK",
        stop.on.error = TRUE,
        progressbar = FALSE
    )
    parallelResult <- do.call(
        permuteRewiringScores,
        c(arguments, list(BPPARAM = snow))
    )

    expect_identical(parallelResult, serialResult)
    expect_false(BiocParallel::bpisup(snow))
})

test_that("permutation progress is reported at the requested interval", {
    candidateEdge <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )

    messages <- capture.output(
        result <- suppressWarnings(permuteRewiringScores(
            analysisData = exampleAnalysisData(),
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            candidateEdgeTable = candidateEdge,
            minimumAbsoluteCorrelation = 0,
            adjustedPValueThreshold = NULL,
            nPermutations = 3L,
            seed = 1L,
            verbose = TRUE,
            progressEvery = 2L
        )),
        type = "message"
    )
    expect_true(is.list(result))
    messages <- paste(messages, collapse = "\n")

    expect_match(messages, "Completed 2 of 3 permutations")
    expect_match(messages, "Completed 3 of 3 permutations")
})

test_that("seed NULL advances the legacy label-permutation stream", {
    analysisData <- exampleAnalysisData()
    sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData))
    permuteLabels <- getFromNamespace(".permuteGroupLabels", "CorNetto")
    candidateEdge <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )

    set.seed(29)
    invisible(lapply(seq_len(3L), function(i) {
        permuteLabels(
            sampleData = sampleData,
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered")
        )
    }))
    expectedState <- .Random.seed

    set.seed(29)
    suppressWarnings(permuteRewiringScores(
        analysisData = analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("PASC", "Recovered"),
        candidateEdgeTable = candidateEdge,
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        nPermutations = 3L,
        seed = NULL
    ))

    expect_identical(.Random.seed, expectedState)
})

test_that("permutation parallel arguments are validated", {
    expect_error(
        suppressWarnings(permuteRewiringScores(
            analysisData = exampleAnalysisData(),
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            nPermutations = 1L,
            BPPARAM = list()
        )),
        "BiocParallelParam"
    )
    expect_error(
        suppressWarnings(permuteRewiringScores(
            analysisData = exampleAnalysisData(),
            groupColumn = "clinicalGroup",
            groupLevels = c("PASC", "Recovered"),
            nPermutations = 1L,
            progressEvery = 0L
        )),
        "`progressEvery`"
    )
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
