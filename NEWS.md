# CorNetto 0.99.0

* Submitted to Bioconductor.

* First release. CorNetto builds knowledge-guided multi-omic correlation
  networks from normalized abundance data held in a `MultiAssayExperiment`,
  and provides:

    * group-specific correlation networks, and pre-analysis summaries of
      sample counts, missingness and feature counts;
    * two-group differential correlation by the Fisher z-difference test,
      either over all within-omic pairs or over candidate edges supplied by a
      prior-knowledge network;
    * node-level rewiring scores, with permutation reference distributions
      that report whether the result is a randomization p-value or a
      conditional ranking;
    * pathway-focused subnetworks, `igraph` conversion, and Cytoscape-ready
      node and edge export.
