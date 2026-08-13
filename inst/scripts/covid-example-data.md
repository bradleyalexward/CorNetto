# Provenance of `inst/extdata/covidData`

This directory documents how the small COVID-19 example files shipped in
`inst/extdata/covidData` were derived. It is a description rather than a
runnable script, because the source data are not public (see Reproducibility
below).

## Source

An unpublished in-house multi-omic COVID-19 cohort collected at UCLouvain and
Cliniques Universitaires Saint-Luc, Brussels. The cohort has no public
repository accession at the time of writing. Participants were sampled at
several visits and classified into severity groups; the full matrices carry
transcriptomic, proteomic and metabolomic abundances that were normalized
before CorNetto ever sees them. CorNetto applies no further transformation.

Contact the package maintainer (see `DESCRIPTION`) for enquiries about the
underlying cohort. An accession will be added here if the cohort is published.

## Shipped files

| File | Contents |
| --- | --- |
| `patientInformation.csv` | 24 rows; columns `Visit_ID`, `Group`, `Visit` |
| `rnaMatrix.csv` | 51 features x 24 samples |
| `proteinMatrix.csv` | 63 features x 24 samples |
| `metaboliteMatrix.csv` | 4 features x 24 samples |

All three matrices carry feature identifiers as row names and sample
identifiers as column names, are complete (no missing values), and hold
normalized log-scale abundances.

## Subsetting

The packaged subset was produced from the full matrices as follows.

1. Restrict `patientInformation.csv` to `Visit == "Visit 1"` and
   `Group %in% c("COVID Moderate", "COVID Severe")`, keeping only participants
   present in all three abundance matrices.
2. Take the first 12 participants of each of the two groups in file order,
   giving 24 samples balanced across severity.
3. Retain a fixed list of features per assay, chosen once during package
   development to keep the example small while leaving enough variable
   features for the vignette's filtering and correlation steps. Selection did
   not depend on any external interaction database, and the retained feature
   identifiers are exactly the row names of the three shipped matrices.
4. Replace participant identifiers with sequential anonymous labels
   `COVID001` to `COVID024`, applied consistently across the patient table and
   all three matrices.
5. Reduce the patient table to `Visit_ID`, `Group` and `Visit`.

## Reproducibility

The full matrices are not distributed with CorNetto: they would exceed the
Bioconductor source-package size limit, and the cohort is not yet public.
Regenerating these files therefore requires access to the source cohort, after
which the five steps above reproduce them exactly. No random numbers are
involved, so the subsetting is deterministic.

## Scope

This subset exists to make the COVID vignette and the package tests runnable
during `R CMD check`. It is far too small for biological inference, and the
candidate edges used alongside it in the vignette are synthetic (see
`make-synthetic-priors.R`). Nothing derived from these files should be read as
a biological result.
