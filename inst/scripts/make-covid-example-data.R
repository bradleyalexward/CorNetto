## This script documents how the small COVID-19 example files in
## inst/extdata/covidData were derived from the local full matrices used
## during CorNetto development.
##
## The full COVID matrices are intentionally not shipped with the package
## because they make the source package too large for Bioconductor. To
## regenerate the small example, set CORNETTO_COVID_FULL_DIR to a directory
## containing:
## - patientInformation.csv
## - rnaMatrix.csv
## - proteinMatrix.csv
## - metaboliteMatrix.csv
##
## The selected feature identifiers are fixed in
## inst/extdata/covidData/featureSelection.csv. Selection does not depend on an
## external interaction database. This keeps regeneration separate from the
## synthetic candidate networks used by the vignette.
##
## Optionally set CORNETTO_COVID_EXAMPLE_DIR to choose the output directory.
## By default, output is written to inst/extdata/covidData.

scriptArg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(scriptArg)) {
    scriptPath <- sub("^--file=", "", scriptArg[[1L]])
} else {
    scriptPath <- file.path("inst", "scripts", "make-covid-example-data.R")
}
packageRoot <- normalizePath(file.path(dirname(scriptPath), "..", ".."))

sourceDir <- Sys.getenv(
    "CORNETTO_COVID_FULL_DIR",
    file.path(packageRoot, "Folder not to be included in final package", "covidDataFull")
)
outputDir <- Sys.getenv(
    "CORNETTO_COVID_EXAMPLE_DIR",
    file.path(packageRoot, "inst", "extdata", "covidData")
)
samplesPerGroup <- as.integer(Sys.getenv("CORNETTO_COVID_SAMPLES_PER_GROUP", "12"))
featureSelectionPath <- file.path(
    packageRoot,
    "inst",
    "extdata",
    "covidData",
    "featureSelection.csv"
)

requiredFiles <- file.path(
    sourceDir,
    c("patientInformation.csv", "rnaMatrix.csv", "proteinMatrix.csv", "metaboliteMatrix.csv")
)
if (!all(file.exists(requiredFiles))) {
    stop(
        "Full COVID source files were not found in `", sourceDir, "`.",
        call. = FALSE
    )
}

readCsvHeader <- function(filePath) {
    headerLine <- readLines(filePath, n = 1L, warn = FALSE)
    scan(
        text = headerLine,
        what = character(),
        sep = ",",
        quote = "\"",
        quiet = TRUE
    )
}

readFeatureIdentifiers <- function(filePath) {
    header <- readCsvHeader(filePath)
    columnClasses <- c("character", rep("NULL", length(header) - 1L))
    featureTable <- read.csv(
        filePath,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = columnClasses
    )
    featureTable[[1L]]
}

readSelectedMatrix <- function(filePath, sampleIds, featureIds) {
    header <- readCsvHeader(filePath)
    columnClasses <- ifelse(header %in% sampleIds, "numeric", "NULL")
    columnClasses[[1L]] <- "character"

    matrixData <- read.csv(
        filePath,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = columnClasses
    )
    names(matrixData)[[1L]] <- "featureIdentifier"
    matrixData <- matrixData[
        matrixData$featureIdentifier %in% featureIds,
        ,
        drop = FALSE
    ]
    rownames(matrixData) <- matrixData$featureIdentifier
    matrixData$featureIdentifier <- NULL

    featureIds <- featureIds[featureIds %in% rownames(matrixData)]
    matrixData[featureIds, sampleIds, drop = FALSE]
}

patientInformation <- read.csv(
    file.path(sourceDir, "patientInformation.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

matrixPaths <- list(
    RNA = file.path(sourceDir, "rnaMatrix.csv"),
    Protein = file.path(sourceDir, "proteinMatrix.csv"),
    Metabolite = file.path(sourceDir, "metaboliteMatrix.csv")
)

matrixSamples <- lapply(matrixPaths, function(path) readCsvHeader(path)[-1L])
availableSamples <- Reduce(intersect, matrixSamples)

covidVisitOne <- patientInformation[
    patientInformation$Visit == "Visit 1" &
        patientInformation$Group %in% c("COVID Moderate", "COVID Severe") &
        patientInformation$Visit_ID %in% availableSamples,
    ,
    drop = FALSE
]

selectedSamples <- unlist(
    lapply(
        c("COVID Moderate", "COVID Severe"),
        function(groupName) {
            head(covidVisitOne$Visit_ID[covidVisitOne$Group == groupName], samplesPerGroup)
        }
    ),
    use.names = FALSE
)

assayFeatures <- lapply(matrixPaths, readFeatureIdentifiers)
featureSelection <- read.csv(
    featureSelectionPath,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
requiredSelectionColumns <- c("assayName", "featureIdentifier")
if (!all(requiredSelectionColumns %in% names(featureSelection))) {
    stop(
        "The feature-selection manifest must contain `assayName` and ",
        "`featureIdentifier`.",
        call. = FALSE
    )
}
if (anyDuplicated(featureSelection) ||
    anyNA(featureSelection$assayName) ||
    anyNA(featureSelection$featureIdentifier) ||
    any(!nzchar(featureSelection$assayName)) ||
    any(!nzchar(featureSelection$featureIdentifier))) {
    stop("The feature-selection manifest contains invalid rows.", call. = FALSE)
}
unknownAssays <- setdiff(unique(featureSelection$assayName), names(assayFeatures))
if (length(unknownAssays)) {
    stop(
        "Unknown assays in the feature-selection manifest: ",
        paste(unknownAssays, collapse = ", "),
        call. = FALSE
    )
}

selectedFeatures <- lapply(names(assayFeatures), function(assayName) {
    selected <- featureSelection$featureIdentifier[
        featureSelection$assayName == assayName
    ]
    missingFeatures <- setdiff(selected, assayFeatures[[assayName]])
    if (length(missingFeatures)) {
        stop(
            "Selected features are absent from the full ", assayName,
            " matrix: ", paste(missingFeatures, collapse = ", "),
            call. = FALSE
        )
    }
    selected
})
names(selectedFeatures) <- names(assayFeatures)

dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

idMap <- setNames(
    sprintf("COVID%03d", seq_along(selectedSamples)),
    selectedSamples
)
patientSubset <- patientInformation[
    match(selectedSamples, patientInformation$Visit_ID),
    c("Visit_ID", "Group", "Visit"),
    drop = FALSE
]
patientSubset$Visit_ID <- unname(idMap[patientSubset$Visit_ID])
write.csv(
    patientSubset,
    file.path(outputDir, "patientInformation.csv"),
    row.names = FALSE,
    na = ""
)

rnaMatrix <- readSelectedMatrix(matrixPaths$RNA, selectedSamples, selectedFeatures$RNA)
proteinMatrix <- readSelectedMatrix(matrixPaths$Protein, selectedSamples, selectedFeatures$Protein)
metaboliteMatrix <- readSelectedMatrix(
    matrixPaths$Metabolite,
    selectedSamples,
    selectedFeatures$Metabolite
)

colnames(rnaMatrix) <- unname(idMap[colnames(rnaMatrix)])
colnames(proteinMatrix) <- unname(idMap[colnames(proteinMatrix)])
colnames(metaboliteMatrix) <- unname(idMap[colnames(metaboliteMatrix)])

write.csv(rnaMatrix, file.path(outputDir, "rnaMatrix.csv"))
write.csv(proteinMatrix, file.path(outputDir, "proteinMatrix.csv"))
write.csv(metaboliteMatrix, file.path(outputDir, "metaboliteMatrix.csv"))

message(
    "Wrote COVID example data with ",
    length(selectedSamples),
    " samples, ",
    nrow(rnaMatrix),
    " RNA features, ",
    nrow(proteinMatrix),
    " protein features, and ",
    nrow(metaboliteMatrix),
    " metabolite features."
)
