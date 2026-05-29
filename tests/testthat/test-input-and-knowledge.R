test_that("example analysis data is a valid MultiAssayExperiment", {
    analysisData <- exampleAnalysisData()
    expect_true(methods::is(analysisData, "MultiAssayExperiment"))

    expect_message(
        validated <- validateAnalysisData(analysisData),
        "structured correctly"
    )
    expect_true(methods::is(validated, "MultiAssayExperiment"))
    expect_silent(validateAnalysisData(analysisData, quiet = TRUE))

    summary <- summarizeAnalysisData(validated, groupColumn = "clinicalGroup", quiet = TRUE)
    expect_true(all(c("assays", "samples", "warnings") %in% names(summary)))
    expect_equal(nrow(summary$assays), 3)
    expect_true(any(grepl("fewer than 10", summary$warnings)))
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

    expect_error(
        filterFeatures(analysisData, assayNames = "protein", topVariableFeatures = -1),
        "positive integer"
    )
    expect_error(
        filterFeatures(analysisData, assayNames = "protein", minimumVariance = -0.1),
        "non-negative"
    )

    assayMatrix <- matrix(
        c(
            1, 1, 1, 1,
            1, 1.1, 1.2, 1.3,
            1, 2, 3, 4
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(
            c("constant", "lowVariance", "highVariance"),
            paste0("S", 1:4)
        )
    )
    sampleData <- data.frame(
        sampleId = paste0("S", 1:4),
        clinicalGroup = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )
    analysisData <- createAnalysisData(list(rna = assayMatrix), sampleData)
    filtered <- filterFeatures(analysisData, removeZeroVariance = TRUE)
    filteredMatrix <- SummarizedExperiment::assay(
        MultiAssayExperiment::experiments(filtered)[["rna"]]
    )
    expect_equal(rownames(filteredMatrix), c("lowVariance", "highVariance"))

    filtered <- filterFeatures(analysisData, minimumVariance = 0.5)
    filteredMatrix <- SummarizedExperiment::assay(
        MultiAssayExperiment::experiments(filtered)[["rna"]]
    )
    expect_equal(rownames(filteredMatrix), "highVariance")
    expect_error(
        filterFeatures(analysisData, removeZeroVariance = NA),
        "removeZeroVariance"
    )
})

test_that("sample filtering preserves assay-specific overlap", {
    analysisData <- exampleAnalysisData()
    experimentList <- MultiAssayExperiment::experiments(analysisData)
    experimentList[["protein"]] <- experimentList[["protein"]][
        ,
        setdiff(colnames(experimentList[["protein"]]), "R1"),
        drop = FALSE
    ]
    MultiAssayExperiment::experiments(analysisData) <- experimentList
    analysisData <- validateAnalysisData(analysisData, quiet = TRUE)

    filtered <- filterSamples(
        analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = "Recovered"
    )

    expect_equal(
        rownames(MultiAssayExperiment::colData(filtered)),
        paste0("R", 1:4)
    )

    filteredExperiments <- MultiAssayExperiment::experiments(filtered)
    proteinAssay <- SummarizedExperiment::assay(
        filteredExperiments[["protein"]]
    )
    transcriptAssay <- SummarizedExperiment::assay(
        filteredExperiments[["transcript"]]
    )
    metaboliteAssay <- SummarizedExperiment::assay(
        filteredExperiments[["metabolite"]]
    )

    expect_equal(ncol(proteinAssay), 3L)
    expect_equal(ncol(transcriptAssay), 4L)
    expect_equal(
        ncol(metaboliteAssay),
        4L
    )
})

test_that("sample filtering validates explicit sample identifiers", {
    analysisData <- exampleAnalysisData()
    filtered <- filterSamples(analysisData, sampleIds = c("R1", "P1"))
    expect_equal(
        rownames(MultiAssayExperiment::colData(filtered)),
        c("R1", "P1")
    )

    expect_error(
        filterSamples(analysisData),
        "Supply `sampleIds`, `groupColumn`, or both"
    )
    expect_error(
        filterSamples(analysisData, sampleIds = "missing"),
        "not found in `colData`"
    )
})

test_that("sampleId metadata fallback and strict numeric assay parsing work", {
    createAnalysisDataArguments <- names(formals(createAnalysisData))
    expect_false("featureIdentifierColumns" %in% createAnalysisDataArguments)
    expect_false("featureNameColumns" %in% createAnalysisDataArguments)

    assayMatrix <- matrix(
        seq_len(12),
        nrow = 3,
        dimnames = list(paste0("F", 1:3), paste0("S", 1:4))
    )
    sampleData <- data.frame(
        sampleId = paste0("S", 1:4),
        clinicalGroup = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    analysisData <- createAnalysisData(list(rna = assayMatrix), sampleData)
    expect_equal(
        rownames(MultiAssayExperiment::colData(analysisData)),
        paste0("S", 1:4)
    )

    assayTable <- data.frame(
        S1 = c(1, 2),
        S2 = c(3, 4),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    expect_error(
        createAnalysisData(list(rna = assayTable), sampleData),
        "row names representing feature identifiers"
    )

    badTable <- data.frame(
        S1 = c("1", "bad"),
        S2 = c("2", "3"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rownames(badTable) <- c("F1", "F2")
    expect_error(
        createAnalysisData(list(rna = badTable), sampleData),
        "Non-numeric values"
    )
})

test_that("example knowledge network validates and can be combined", {
    standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
    getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
    knowledgeNetwork <- exampleKnowledgeNetwork()
    expect_true(methods::is(knowledgeNetwork, "DataFrame"))
    expect_true(all(standardColumns %in% names(knowledgeNetwork)))
    expect_gt(
        nrow(getNodeTable(knowledgeNetwork, fallbackToEdges = FALSE)),
        0
    )
    expect_message(
        validatedKnowledgeNetwork <- validateKnowledgeNetwork(knowledgeNetwork),
        "structured correctly"
    )
    expect_true(methods::is(validatedKnowledgeNetwork, "DataFrame"))
    expect_silent(validateKnowledgeNetwork(knowledgeNetwork, quiet = TRUE))

    analysisData <- exampleAnalysisData()
    analysisData <- validateKnowledgeNetwork(
        knowledgeNetwork,
        analysisData = analysisData,
        resultName = "storedKnowledgeNetwork",
        storeResult = TRUE,
        quiet = TRUE
    )
    expect_identical(
        getCorNettoResult(
            analysisData,
            "knowledgeNetworks",
            "storedKnowledgeNetwork"
        ),
        validatedKnowledgeNetwork
    )

    combined <- combineKnowledgeNetworks(
        knowledgeNetwork,
        knowledgeNetwork,
        removeDuplicates = TRUE
    )
    expect_true(methods::is(combined, "DataFrame"))
    expect_equal(nrow(combined), nrow(knowledgeNetwork))
    expect_gt(nrow(getNodeTable(combined, fallbackToEdges = FALSE)), 0)

    networkSummary <- summarizeKnowledgeNetwork(
        combined,
        analysisData = exampleAnalysisData(),
        quiet = TRUE
    )
    expect_true(all(c("overall", "assayPairs", "warnings") %in% names(networkSummary)))
    expect_equal(networkSummary$overall$totalEdges, nrow(combined))
    expect_gt(networkSummary$overall$measuredAssayEdges, 0)
})
