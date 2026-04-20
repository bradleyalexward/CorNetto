## This script documents the synthetic example assets shipped with
## CorNetto. The current package build stores these example resources as
## plain text files in inst/extdata so they remain inspectable and do not
## require local R package build tooling to regenerate.
##
## The synthetic data are intentionally small and deterministic:
## - 8 samples split into Recovered and PASC groups
## - 3 normalized assays: protein, transcript, metabolite
## - a compact knowledge network spanning within-omic and cross-omic
##   priors plus a drug-target edge
## - a small seed-node set for pathway-focused workflow examples
