## Build the synthetic candidate-edge files used by the COVID vignette.
##
## The edges are generated from feature identifiers in the packaged matrices.
## They are arbitrary software examples, not biological assertions and not
## extracts from external interaction databases. No random numbers are used.

scriptArg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
scriptPath <- if (length(scriptArg)) {
    sub("^--file=", "", scriptArg[[1L]])
} else {
    file.path("inst", "scripts", "make-synthetic-priors.R")
}
packageRoot <- normalizePath(file.path(dirname(scriptPath), "..", ".."))
covidDir <- file.path(packageRoot, "inst", "extdata", "covidData")
outputDir <- file.path(packageRoot, "inst", "extdata", "priorNetworks")

patientData <- read.csv(
    file.path(covidDir, "patientInformation.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
selectedSamples <- patientData$Visit_ID[
    patientData$Group %in% c("COVID Moderate", "COVID Severe") &
        patientData$Visit == "Visit 1"
]

retainedFeatures <- function(fileName) {
    abundance <- read.csv(
        file.path(covidDir, fileName),
        row.names = 1L,
        check.names = FALSE
    )
    sampleIds <- selectedSamples[selectedSamples %in% colnames(abundance)]
    featureVariance <- apply(
        abundance[, sampleIds, drop = FALSE],
        1L,
        stats::var,
        na.rm = TRUE
    )
    names(featureVariance)[
        is.finite(featureVariance) & featureVariance >= 0.01
    ]
}

features <- list(
    RNA = retainedFeatures("rnaMatrix.csv"),
    Protein = retainedFeatures("proteinMatrix.csv"),
    Metabolite = retainedFeatures("metaboliteMatrix.csv")
)
minimumFeatureCounts <- c(RNA = 8L, Protein = 12L, Metabolite = 4L)
if (any(lengths(features) < minimumFeatureCounts)) {
    stop(
        "Synthetic-prior generation requires at least 8 RNA, 12 protein, ",
        "and 4 metabolite features after variance filtering.",
        call. = FALSE
    )
}

makeEdges <- function(
    from,
    to,
    fromAssay,
    toAssay,
    edgeType,
    edgeDirection,
    weightBase,
    weightStep,
    weightCycle
) {
    n <- length(from)
    stopifnot(length(to) == n)
    data.frame(
        from = from,
        to = to,
        fromName = from,
        toName = to,
        edgeWeight = round(
            weightBase + ((seq_len(n) - 1L) %% weightCycle) * weightStep,
            3L
        ),
        edgeDirection = edgeDirection,
        edgeType = edgeType,
        fromOmic = fromAssay,
        toOmic = toAssay,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

withinAssay <- do.call(
    rbind,
    lapply(names(features), function(assayName) {
        featureIds <- features[[assayName]]
        nextFeature <- featureIds[c(seq.int(2L, length(featureIds)), 1L)]
        makeEdges(
            featureIds,
            nextFeature,
            assayName,
            assayName,
            "syntheticWithinAssay",
            "undirected",
            0.55,
            0.05,
            9L
        )
    })
)
rownames(withinAssay) <- NULL

rna <- features$RNA
protein <- features$Protein
metabolite <- features$Metabolite

rnaProtein <- makeEdges(
    rna,
    protein[((seq_along(rna) - 1L) %% length(protein)) + 1L],
    "RNA",
    "Protein",
    "syntheticCrossAssay",
    "from_to",
    0.60,
    0.04,
    8L
)
proteinMetabolite <- makeEdges(
    protein[seq_len(12L)],
    metabolite[((seq_len(12L) - 1L) %% length(metabolite)) + 1L],
    "Protein",
    "Metabolite",
    "syntheticCrossAssay",
    "from_to",
    0.62,
    0.04,
    8L
)
rnaMetabolite <- makeEdges(
    rna[seq_len(8L)],
    metabolite[((seq_len(8L) - 1L) %% length(metabolite)) + 1L],
    "RNA",
    "Metabolite",
    "syntheticCrossAssay",
    "from_to",
    0.64,
    0.04,
    8L
)
crossAssay <- rbind(rnaProtein, proteinMetabolite, rnaMetabolite)
rownames(crossAssay) <- NULL

decoyEdges <- data.frame(
    from = c(
        rna[[1L]], "SYN_RNA_01", protein[[1L]], "SYN_PROTEIN_01",
        metabolite[[1L]], "SYN_METABOLITE_01", "SYN_CONTEXT_01",
        "SYN_CONTEXT_03"
    ),
    to = c(
        "SYN_RNA_01", rna[[2L]], "SYN_PROTEIN_01", protein[[2L]],
        "SYN_METABOLITE_01", "SYN_METABOLITE_02", rna[[3L]],
        "SYN_CONTEXT_04"
    ),
    fromName = c(
        rna[[1L]], "unmeasured RNA 1", protein[[1L]],
        "unmeasured protein 1", metabolite[[1L]],
        "unmeasured metabolite 1", "unmeasured context 1",
        "unmeasured context 3"
    ),
    toName = c(
        "unmeasured RNA 1", rna[[2L]], "unmeasured protein 1",
        protein[[2L]], "unmeasured metabolite 1",
        "unmeasured metabolite 2", rna[[3L]], "unmeasured context 4"
    ),
    edgeWeight = seq(0.50, 0.78, length.out = 8L),
    edgeDirection = rep(c("undirected", "from_to"), 4L),
    edgeType = "syntheticDecoy",
    fromOmic = c(
        "RNA", "RNA", "Protein", "Protein", "Metabolite", "Metabolite",
        "SyntheticContext", "SyntheticContext"
    ),
    toOmic = c(
        "RNA", "RNA", "Protein", "Protein", "Metabolite", "Metabolite",
        "RNA", "SyntheticContext"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

allEdges <- rbind(withinAssay, crossAssay, decoyEdges)
measuredKeys <- unlist(
    Map(function(assayName, featureIds) {
        paste(assayName, featureIds, sep = "::")
    }, names(features), features),
    use.names = FALSE
)
measuredEndpoints <-
    (paste(allEdges$fromOmic, allEdges$from, sep = "::") %in% measuredKeys) +
    (paste(allEdges$toOmic, allEdges$to, sep = "::") %in% measuredKeys)
stopifnot(
    any(measuredEndpoints == 2L),
    any(measuredEndpoints == 1L),
    any(measuredEndpoints == 0L),
    length(unique(allEdges$edgeType)) == 3L
)

dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
outputs <- list(
    synthetic_within_assay.csv = withinAssay,
    synthetic_cross_assay.csv = crossAssay,
    synthetic_decoy_edges.csv = decoyEdges
)
for (fileName in names(outputs)) {
    filePath <- file.path(outputDir, fileName)
    write.csv(
        outputs[[fileName]],
        filePath,
        row.names = FALSE,
        quote = FALSE,
        fileEncoding = "UTF-8"
    )
    message(fileName, ": ", nrow(outputs[[fileName]]), " edges")
}
