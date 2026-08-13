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
    expect_equal(summary$samples$overall$sampleCount, 20L)
    expect_equal(summary$samples$groups$sampleCount, c(10L, 10L))
    # The example data sits at the default threshold, so nothing should fire.
    expect_false(any(grepl("fewer than", summary$warnings)))

    strictSummary <- summarizeAnalysisData(
        validated,
        groupColumn = "clinicalGroup",
        lowSampleWarningThreshold = 15L,
        quiet = TRUE
    )
    expect_true(any(grepl("fewer than 15", strictSummary$warnings)))
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
        "finite integer"
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

test_that("feature filtering handles assay-specific named-list arguments", {
    analysisData <- exampleAnalysisData()
    experimentList <- MultiAssayExperiment::experiments(analysisData)
    originalFeatureCounts <- vapply(
        experimentList,
        function(assayObject) nrow(SummarizedExperiment::assay(assayObject)),
        integer(1L)
    )
    proteinFeatures <- rownames(SummarizedExperiment::assay(
        experimentList[["protein"]]
    ))[seq_len(2L)]

    filteredSubset <- filterFeatures(
        analysisData,
        featureSubset = list(protein = proteinFeatures)
    )
    subsetFeatureCounts <- vapply(
        MultiAssayExperiment::experiments(filteredSubset),
        function(assayObject) nrow(SummarizedExperiment::assay(assayObject)),
        integer(1L)
    )
    expect_equal(subsetFeatureCounts[["protein"]], 2L)
    expect_equal(
        subsetFeatureCounts[c("transcript", "metabolite")],
        originalFeatureCounts[c("transcript", "metabolite")]
    )

    filteredTop <- filterFeatures(
        analysisData,
        topVariableFeatures = list(protein = 3L)
    )
    topFeatureCounts <- vapply(
        MultiAssayExperiment::experiments(filteredTop),
        function(assayObject) nrow(SummarizedExperiment::assay(assayObject)),
        integer(1L)
    )
    expect_equal(topFeatureCounts[["protein"]], 3L)
    expect_equal(
        topFeatureCounts[c("transcript", "metabolite")],
        originalFeatureCounts[c("transcript", "metabolite")]
    )

    expect_silent(filterFeatures(
        analysisData,
        minimumVariance = list(protein = 0)
    ))
    expect_error(
        filterFeatures(analysisData, featureSubset = list(proteinFeatures)),
        "named lists"
    )
})

test_that("feature filtering drops assays with no retained features", {
    analysisData <- exampleAnalysisData()
    filtered <- NULL
    expect_warning(
        filtered <- filterFeatures(
            analysisData = analysisData,
            assayNames = c("protein", "transcript", "metabolite"),
            featureSubset = list(
                protein = c("P1", "P2"),
                transcript = c("T1", "T2"),
                metabolite = "not_measured"
            )
        ),
        "dropping this assay"
    )

    filteredAssays <- names(MultiAssayExperiment::experiments(filtered))
    expect_true(all(c("protein", "transcript") %in% filteredAssays))
    expect_false("metabolite" %in% filteredAssays)
    expect_equal(
        nrow(SummarizedExperiment::assay(
            MultiAssayExperiment::experiments(filtered)[["protein"]]
        )),
        2L
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
        paste0("R", 1:10)
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

    # protein lost R1 above, the other assays keep all ten Recovered samples
    expect_equal(ncol(proteinAssay), 9L)
    expect_equal(ncol(transcriptAssay), 10L)
    expect_equal(ncol(metaboliteAssay), 10L)
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

test_that("packaged synthetic priors cover measured and decoy edges", {
    priorDir <- system.file("extdata", "priorNetworks", package = "CorNetto")
    fileNames <- c(
        "synthetic_within_assay.csv",
        "synthetic_cross_assay.csv",
        "synthetic_decoy_edges.csv"
    )
    expect_true(nzchar(priorDir))
    expect_true(all(file.exists(file.path(priorDir, fileNames))))

    columnMap <- c(
        fromFeatureIdentifier = "from",
        toFeatureIdentifier = "to",
        fromFeatureName = "fromName",
        toFeatureName = "toName",
        fromAssayName = "fromOmic",
        toAssayName = "toOmic",
        edgeType = "edgeType",
        edgeDirection = "edgeDirection",
        evidenceScore = "edgeWeight"
    )
    networks <- lapply(fileNames, function(fileName) {
        readKnowledgeNetwork(
            file.path(priorDir, fileName),
            columnMapping = columnMap,
            knowledgeSource = paste("synthetic", fileName)
        )
    })
    combined <- do.call(
        combineKnowledgeNetworks,
        c(networks, list(removeDuplicates = TRUE))
    )

    expect_setequal(
        unique(combined$edgeType),
        c("syntheticWithinAssay", "syntheticCrossAssay", "syntheticDecoy")
    )
    expect_true(any(combined$isDirected))
    expect_true(any(!combined$isDirected))

    covidDir <- system.file("extdata", "covidData", package = "CorNetto")
    assayFiles <- c(
        RNA = "rnaMatrix.csv",
        Protein = "proteinMatrix.csv",
        Metabolite = "metaboliteMatrix.csv"
    )
    measuredKeys <- unlist(
        lapply(names(assayFiles), function(assayName) {
            abundance <- read.csv(
                file.path(covidDir, assayFiles[[assayName]]),
                row.names = 1L,
                check.names = FALSE
            )
            paste(assayName, rownames(abundance), sep = "::")
        }),
        use.names = FALSE
    )
    measured <- filterNetworkByNodes(
        combined,
        nodes = measuredKeys,
        mode = "both",
        nodeIdType = "key"
    )

    expect_lt(nrow(measured), nrow(combined))
    expect_true(all(
        paste(measured$fromAssayName, measured$fromFeatureIdentifier, sep = "::") %in%
            measuredKeys
    ))
    expect_true(all(
        paste(measured$toAssayName, measured$toFeatureIdentifier, sep = "::") %in%
            measuredKeys
    ))
})

test_that("analysis input requires explicit sample and assay identifiers", {
    assay <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(c("F1", "F2"), paste0("S", 1:4))
    )
    automaticRows <- data.frame(group = rep(c("A", "B"), each = 2))
    expect_error(
        createAnalysisData(list(protein = assay), automaticRows),
        "explicit row names"
    )

    sampleData <- data.frame(
        sampleId = paste0("S", 1:4),
        group = rep(c("A", "B"), each = 2)
    )
    duplicateNames <- list(assay, assay)
    names(duplicateNames) <- c("protein", "protein")
    expect_error(
        createAnalysisData(duplicateNames, sampleData),
        "unique"
    )
})

test_that("analysis input rejects NaN and complex abundances", {
    sampleData <- data.frame(sampleId = c("S1", "S2"), group = c("A", "B"))
    nanAssay <- matrix(
        c(1, NaN, 3, 4),
        nrow = 2,
        dimnames = list(c("F1", "F2"), c("S1", "S2"))
    )
    complexAssay <- matrix(
        c(1 + 1i, 2 + 0i, 3 + 0i, 4 + 0i),
        nrow = 2,
        dimnames = list(c("F1", "F2"), c("S1", "S2"))
    )

    expect_error(createAnalysisData(list(protein = nanAssay), sampleData), "Non-finite")
    expect_error(createAnalysisData(list(protein = complexAssay), sampleData), "Complex")
})

test_that("standard edge coercion preserves numeric factor values", {
    edge <- data.frame(
        fromFeatureIdentifier = "A",
        toFeatureIdentifier = "B",
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeWeight = factor("-2"),
        pValue = factor("0.025"),
        stringsAsFactors = TRUE
    )

    result <- coerceEdgeTable(edge)
    expect_equal(result$edgeWeight, -2)
    expect_equal(result$pValue, 0.025)
})

test_that("standard edge coercion preserves logical factor values", {
    edge <- data.frame(
        fromFeatureIdentifier = c("A", "B"),
        toFeatureIdentifier = c("B", "C"),
        fromAssayName = "protein",
        toAssayName = "protein",
        isDirected = factor(c("TRUE", "FALSE"))
    )

    result <- coerceEdgeTable(edge)
    expect_identical(result$isDirected, c(TRUE, FALSE))
})

test_that("knowledge networks reject blank endpoints", {
    network <- data.frame(
        fromFeatureIdentifier = " ",
        toFeatureIdentifier = "B",
        fromAssayName = "protein",
        toAssayName = "protein"
    )
    expect_error(validateKnowledgeNetwork(network, quiet = TRUE), "blank")
})
