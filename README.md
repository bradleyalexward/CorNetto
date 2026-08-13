# CorNetto <img src="man/figures/logo.png" align="right" width="130" alt="" />

CorNetto is a Bioconductor-style R package for knowledge-guided
multi-omic correlation network analysis. It is designed for normalized
transcriptomic, proteomic, and metabolomic abundance data stored in a
`MultiAssayExperiment`.

The first package version focuses on:

- group-specific correlation networks
- pre-analysis summaries of sample counts, missingness, and feature counts
- two-group differential correlation testing, over all pairs in an assay or
  over candidate edges taken from a prior-knowledge network
- rewiring scores for perturbed nodes
- permutation-based assessment of rewiring scores
- pathway-focused subnetworks
- graph construction and Cytoscape-ready export

## Workflow

Green boxes are functions, blue boxes are the objects they return, and orange
boxes are inputs and the `MultiAssayExperiment` that carries results between
steps.

![CorNetto workflow](man/figures/flowchart.png)

## Installation

CorNetto is under review for Bioconductor. Once accepted, install it with:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CorNetto")
```

Until then, the development version can be installed from GitHub:

```r
BiocManager::install("bradleyalexward/CorNetto")
```

## Getting started

```r
library(CorNetto)

analysisData <- exampleAnalysisData()
knowledgeNetwork <- exampleKnowledgeNetwork()

# With assayName omitted, every assay is tested densely within-omic.
# Protein-protein and transcript-transcript pairs are included, for example,
# but protein-transcript pairs are not generated.
differentialResults <- testDifferentialCorrelation(
    analysisData,
    groupColumn = "clinicalGroup",
    groupLevels = c("Recovered", "PASC"),
    storeResult = FALSE
)

differentialNetwork <- createDifferentialCorrelationNetwork(differentialResults)
rewiringScores <- calculateRewiringScores(differentialNetwork)
```

Long permutation runs can use a BiocParallel backend. `SnowParam` is the
portable choice on Windows:

```r
backend <- BiocParallel::SnowParam(
    workers = 2,
    type = "SOCK",
    stop.on.error = TRUE,
    progressbar = FALSE
)

rewiringValidation <- permuteRewiringScores(
    analysisData,
    candidateEdgeTable = knowledgeNetwork,
    groupColumn = "clinicalGroup",
    groupLevels = c("Recovered", "PASC"),
    minimumAbsoluteCorrelation = 0,
    adjustedPValueThreshold = NULL,
    nPermutations = 999,
    seed = 1,
    verbose = TRUE,
    progressEvery = 300,
    BPPARAM = backend
)
```

See `vignette("CorNetto")` for the full workflow.

## Citation

If you use or build upon CorNetto, please cite it. The full author list,
affiliations and a formatted reference are returned by:

```r
citation("CorNetto")
```

A machine-readable citation is also available in `CITATION.cff`.

## License

CorNetto is released under the Artistic License 2.0.
