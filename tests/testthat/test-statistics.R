## These compare CorNetto's numbers against stats::cor.test, against the
## published formulas, and against hand-worked values. They are the tests that
## would actually catch a wrong statistic, as opposed to a wrong column name.

test_that("dense pearson correlations and p-values match stats::cor.test", {
    set.seed(42)
    m <- matrix(
        rnorm(60),
        nrow = 4,
        dimnames = list(paste0("F", 1:4), paste0("S", 1:15))
    )
    result <- computePearson(m)

    for (i in seq_len(3L)) {
        for (j in seq.int(i + 1L, 4L)) {
            reference <- stats::cor.test(m[i, ], m[j, ], method = "pearson")
            expect_equal(result$correlationMatrix[i, j], unname(reference$estimate))
            expect_equal(result$pValueMatrix[i, j], reference$p.value)
            expect_equal(result$sampleCountMatrix[i, j], 15)
        }
    }
})

test_that("dense pearson p-values match cor.test under pairwise missingness", {
    set.seed(7)
    m <- matrix(
        rnorm(60),
        nrow = 4,
        dimnames = list(paste0("F", 1:4), paste0("S", 1:15))
    )
    m[1, 1:5] <- NA

    result <- computePearson(m)
    keep <- stats::complete.cases(m[1, ], m[2, ])
    reference <- stats::cor.test(m[1, keep], m[2, keep], method = "pearson")

    expect_equal(result$correlationMatrix[1, 2], unname(reference$estimate))
    expect_equal(result$pValueMatrix[1, 2], reference$p.value)
    expect_equal(result$sampleCountMatrix[1, 2], 10)
})

test_that("perfect pearson correlations match cor.test", {
    m <- rbind(
        increasing = 1:6,
        same = 2 * (1:6),
        decreasing = -3 * (1:6)
    )

    result <- computePearson(m)
    positive <- stats::cor.test(m[1, ], m[2, ], method = "pearson")
    negative <- stats::cor.test(m[1, ], m[3, ], method = "pearson")

    expect_equal(result$pValueMatrix[1, 2], positive$p.value)
    expect_equal(result$pValueMatrix[1, 3], negative$p.value)
    expect_equal(result$pValueMatrix[1, 2], 0)
    expect_equal(result$pValueMatrix[1, 3], 0)
})

test_that("spearman correlations match cor.test with the normal approximation", {
    set.seed(11)
    m <- matrix(
        rnorm(45),
        nrow = 3,
        dimnames = list(paste0("F", 1:3), paste0("S", 1:15))
    )
    result <- computePairwise(m, correlationMethod = "spearman")
    reference <- suppressWarnings(
        stats::cor.test(m[1, ], m[2, ], method = "spearman", exact = FALSE)
    )

    expect_equal(result$correlationMatrix[1, 2], unname(reference$estimate))
    expect_equal(result$pValueMatrix[1, 2], reference$p.value)
})

test_that("kendall correlations match cor.test with the normal approximation", {
    set.seed(13)
    m <- matrix(
        rnorm(45),
        nrow = 3,
        dimnames = list(paste0("F", 1:3), paste0("S", 1:15))
    )
    result <- computePairwise(m, correlationMethod = "kendall")
    reference <- suppressWarnings(
        stats::cor.test(m[1, ], m[2, ], method = "kendall", exact = FALSE)
    )

    expect_equal(result$correlationMatrix[1, 2], unname(reference$estimate))
    expect_equal(result$pValueMatrix[1, 2], reference$p.value)
})

test_that("the z-score difference is the Fisher z difference from the literature", {
    # Fisher (1921). Group 1 is subtracted from group 2, so a positive
    # z means the correlation is stronger in group 2. See the sign-convention
    # test below, which is what pins this down.
    expected <- (atanh(0.2) - atanh(0.8)) / sqrt(1 / (20 - 3) + 1 / (20 - 3))
    observed <- fisherDifference(0.8, 20, 0.2, 20, "pearson")

    expect_equal(unname(observed), expected, tolerance = 1e-8)
    expect_equal(
        2 * stats::pnorm(abs(unname(observed)), lower.tail = FALSE),
        2 * stats::pnorm(abs(expected), lower.tail = FALSE)
    )
})

test_that("zScoreDifference is positive when group 2 is more correlated", {
    # The sign convention is load-bearing: edgeWeightMethod = "signedZScore"
    # exposes it directly as the edge weight, while the edgeDirection labels
    # ("gainInGroup1") are group-1 centric. Pin the direction here.
    strongerInGroup2 <- fisherDifference(0.1, 20, 0.9, 20, "pearson")
    strongerInGroup1 <- fisherDifference(0.9, 20, 0.1, 20, "pearson")

    expect_gt(unname(strongerInGroup2), 0)
    expect_lt(unname(strongerInGroup1), 0)
    expect_equal(unname(strongerInGroup2), -unname(strongerInGroup1))
})

test_that("the spearman z difference carries the Fieller variance inflation", {
    # Fieller, Hartley and Pearson (1957) doi:10.1093/biomet/44.3-4.470
    expected <- (atanh(0.2) - atanh(0.8)) / sqrt(1.06 / (20 - 3) + 1.06 / (20 - 3))
    observed <- fisherDifference(0.8, 20, 0.2, 20, "spearman")

    expect_equal(unname(observed), expected, tolerance = 1e-8)
})

test_that("unequal group sizes enter the z difference asymmetrically", {
    expected <- (atanh(0.1) - atanh(0.9)) / sqrt(1 / (12 - 3) + 1 / (30 - 3))
    observed <- fisherDifference(0.9, 12, 0.1, 30, "pearson")

    expect_equal(unname(observed), expected, tolerance = 1e-8)
})

test_that("candidate-edge z-scores and p-values match the hand formula", {
    analysisData <- exampleAnalysisData()
    proteinMatrix <- SummarizedExperiment::assay(
        MultiAssayExperiment::experiments(analysisData)[["protein"]]
    )
    sampleData <- as.data.frame(MultiAssayExperiment::colData(analysisData))
    recovered <- rownames(sampleData)[sampleData$clinicalGroup == "Recovered"]
    pasc <- rownames(sampleData)[sampleData$clinicalGroup == "PASC"]

    candidate <- S4Vectors::DataFrame(
        fromFeatureIdentifier = "P1",
        toFeatureIdentifier = "P2",
        fromAssayName = "protein",
        toAssayName = "protein",
        check.names = FALSE
    )
    result <- testDifferentialCorrelation(
        analysisData,
        groupColumn = "clinicalGroup",
        groupLevels = c("Recovered", "PASC"),
        candidateEdgeTable = candidate,
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        storeResult = FALSE
    )
    expect_equal(nrow(result), 1L)

    r1 <- stats::cor(proteinMatrix["P1", recovered], proteinMatrix["P2", recovered])
    r2 <- stats::cor(proteinMatrix["P1", pasc], proteinMatrix["P2", pasc])
    # Group 2 minus group 1.
    expectedZ <- (atanh(r2) - atanh(r1)) /
        sqrt(1 / (length(recovered) - 3) + 1 / (length(pasc) - 3))

    expect_equal(result$group1CorrelationValue, r1)
    expect_equal(result$group2CorrelationValue, r2)
    expect_equal(result$zScoreDifference, expectedZ, tolerance = 1e-6)
    # P1-P2 is the designed rewiring pair: tightly coupled in Recovered
    # (group 1), decoupled in PASC (group 2), so the z is negative.
    expect_lt(result$zScoreDifference, 0)
    expect_equal(
        result$pValue,
        2 * stats::pnorm(abs(expectedZ), lower.tail = FALSE),
        tolerance = 1e-6
    )
    expect_equal(result$group1SampleCount, length(recovered))
    expect_equal(result$group2SampleCount, length(pasc))
})

test_that("dense differential correlation uses pairwise sample counts", {
    set.seed(29)
    sampleIds <- paste0("S", seq_len(20))
    assay <- matrix(
        rnorm(40),
        nrow = 2,
        dimnames = list(c("F1", "F2"), sampleIds)
    )
    assay[1, 1:4] <- NA_real_
    sampleData <- data.frame(
        sampleId = sampleIds,
        group = rep(c("A", "B"), each = 10),
        stringsAsFactors = FALSE
    )
    analysisData <- createAnalysisData(list(protein = assay), sampleData)

    result <- testDifferentialCorrelation(
        analysisData,
        groupColumn = "group",
        groupLevels = c("A", "B"),
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = NULL,
        storeResult = FALSE
    )

    completeA <- stats::complete.cases(assay[1, 1:10], assay[2, 1:10])
    r1 <- stats::cor(assay[1, 1:10][completeA], assay[2, 1:10][completeA])
    r2 <- stats::cor(assay[1, 11:20], assay[2, 11:20])
    expected <- (atanh(r2) - atanh(r1)) /
        sqrt(1 / (sum(completeA) - 3) + 1 / (10 - 3))

    expect_equal(result$group1SampleCount, 6)
    expect_equal(result$group2SampleCount, 10)
    expect_equal(result$zScoreDifference, expected)
})

test_that("permutation tail probabilities follow Phipson and Smyth", {
    # (1 + r) / (1 + B); doi:10.2202/1544-6115.1585

    nullValues <- c(1, 2, 3, 4, 5)
    expect_equal(tailProbability(3, nullValues, "greater"), (1 + 3) / (1 + 5))
    expect_equal(tailProbability(6, nullValues, "greater"), (1 + 0) / (1 + 5))
    expect_equal(tailProbability(0, nullValues, "greater"), (1 + 5) / (1 + 5))

    # Unscored permutations drop out of the denominator rather than counting as
    # a null draw of zero.
    expect_equal(
        tailProbability(3, c(nullValues, NA, NA), "greater"),
        (1 + 3) / (1 + 5)
    )
    expect_true(is.na(tailProbability(3, c(NA, NA), "greater")))
    expect_true(is.na(tailProbability(NA, nullValues, "greater")))
    expect_error(tailProbability(3, nullValues, "sideways"), "Unknown")
})

test_that("rewiring scores are the incident-edge L2 norm", {
    # A star: N1 carries three edges with z = 3, 4, 12.
    network <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("N1", "N1", "N1"),
        toFeatureIdentifier = c("N2", "N3", "N4"),
        fromFeatureName = c("N1", "N1", "N1"),
        toFeatureName = c("N2", "N3", "N4"),
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeWeight = c(3, 4, 12),
        isDirected = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    scores <- calculateRewiringScores(network, storeResult = FALSE)
    hub <- as.data.frame(scores)[scores$nodeIdentifier == "N1", ]

    expect_equal(hub$totalConnections, 3L)
    expect_equal(hub$rawRewiringScore, sqrt(3^2 + 4^2 + 12^2))
    expect_equal(hub$rawRewiringScore, 13)
    expect_equal(hub$rootMeanSquareRewiringScore, 13 / sqrt(3))

    leaf <- as.data.frame(scores)[scores$nodeIdentifier == "N2", ]
    expect_equal(leaf$totalConnections, 1L)
    expect_equal(leaf$rawRewiringScore, 3)
    expect_equal(leaf$rootMeanSquareRewiringScore, 3)
})

test_that("rewiring scores reject missing and non-finite weights", {
    baseEdge <- data.frame(
        fromFeatureIdentifier = "A",
        toFeatureIdentifier = "B",
        fromAssayName = "protein",
        toAssayName = "protein",
        edgeWeight = NA_real_,
        stringsAsFactors = FALSE
    )

    expect_error(
        calculateRewiringScores(coerceEdgeTable(baseEdge)),
        "finite, non-missing"
    )
    baseEdge$edgeWeight <- Inf
    # A raw table reaches the shared numeric-column validator before the
    # rewiring-specific missing-weight check above.
    expect_error(
        calculateRewiringScores(baseEdge),
        "`edgeWeight` cannot contain non-finite values"
    )
})

test_that("degree-matched z-scores are NA for bins that cannot support them", {
    # Six nodes at the same degree: a real spread, so a real z.
    scored <- computeDegreeMatched(
        rawRewiringScore = c(1, 2, 3, 4, 5, 6),
        totalConnections = rep(3L, 6L)
    )
    expect_false(any(is.na(scored$degreeMatchedZScore)))
    expect_equal(mean(scored$degreeMatchedZScore), 0)
    expect_equal(stats::sd(scored$degreeMatchedZScore), 1)

    # Two nodes in a bin would give +/-0.707 whatever the data says.
    tiny <- computeDegreeMatched(
        rawRewiringScore = c(1, 9),
        totalConnections = c(3L, 3L)
    )
    expect_true(all(is.na(tiny$degreeMatchedZScore)))

    # No spread is not estimable either.
    flat <- computeDegreeMatched(
        rawRewiringScore = rep(5, 6L),
        totalConnections = rep(3L, 6L)
    )
    expect_true(all(is.na(flat$degreeMatchedZScore)))
})

test_that("pairs with fewer than three shared samples do not become NA edges", {
    # Regression test: cor() returns +/-1 through two points while the p-value
    # is NA, so keepIndex used to carry NA and inject all-missing rows.
    set.seed(3)
    m <- matrix(
        rnorm(30),
        nrow = 3,
        dimnames = list(paste0("F", 1:3), paste0("S", 1:10))
    )
    m[1, 3:10] <- NA

    sampleData <- data.frame(
        sampleId = paste0("S", 1:10),
        clinicalGroup = "A",
        stringsAsFactors = FALSE
    )
    analysisData <- createAnalysisData(list(protein = m), sampleData)

    result <- suppressWarnings(createCorrelationNetwork(
        analysisData,
        assayName = "protein",
        minimumAbsoluteCorrelation = 0,
        adjustedPValueThreshold = 1,
        storeResult = FALSE
    ))

    expect_false(any(is.na(result$fromFeatureIdentifier)))
    expect_false(any(is.na(result$toFeatureIdentifier)))
    # F2-F3 is testable; both pairs involving F1 are not.
    expect_equal(nrow(result), 1L)
    expect_equal(sort(c(result$fromFeatureIdentifier, result$toFeatureIdentifier)),
                 c("F2", "F3"))
})

test_that("a feature constant within one group does not crash classification", {
    differentialResults <- coerceEdgeTable(data.frame(
        fromFeatureIdentifier = c("A", "B"),
        toFeatureIdentifier = c("C", "D"),
        fromFeatureName = c("A", "B"),
        toFeatureName = c("C", "D"),
        fromAssayName = "protein",
        toAssayName = "protein",
        group1CorrelationValue = c(0.9, 0.9),
        group2CorrelationValue = c(NA_real_, 0.1),
        zScoreDifference = c(2, 3),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))

    network <- createDifferentialCorrelationNetwork(
        differentialResults,
        minimumAbsoluteCorrelation = 0.3
    )
    expect_equal(nrow(network), 2L)
    expect_true(is.na(network$edgeDirection[[1]]))
    expect_equal(network$edgeDirection[[2]], "gainInGroup1")
})
