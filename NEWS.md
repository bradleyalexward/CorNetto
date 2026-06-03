CHANGES IN VERSION 0.99.1
-------------------------

NEW FEATURES

    o Added summarizeKnowledgeNetwork() for prior-network diagnostics.

SIGNIFICANT USER-VISIBLE CHANGES

    o Simplified the public network API: combineNetworks() now handles
      network integration, getCorNettoResult() is the main result accessor,
      correlationResults() replaces correlation-network-specific accessors,
      and permuteRewiringScores() replaces permuteDifferentialRewiring().

    o Removed readAssayData(), mapFeatureIdentifiers(),
      validateSparseCorrelationEdges(), createCorrelationNetworks(),
      combineCorrelationNetworks(), and createIntegratedNetwork().

    o Replaced the full COVID-19 demonstration matrices with a small
      reproducible example subset suitable for Bioconductor package checks.

    o Updated the differential-correlation to rewiring workflow so
      testDifferentialCorrelation() first applies the absolute-correlation
      prefilter, adjusts p-values on the retained edge universe, and keeps
      edges only when the same group is both sufficiently correlated and
      significant. createDifferentialCorrelationNetwork() no longer applies
      a differential-correlation p-value filter by default and now uses
      signed z-score differences as the default edge weights.

BUG FIXES

    o Fixed storage of named results in MultiAssayExperiment metadata.

CHANGES IN VERSION 0.99.0
-------------------------

NEW FEATURES

    o Initial Bioconductor-oriented package scaffold.

    o Added support for normalized multi-omic inputs through
      MultiAssayExperiment.

    o Implemented knowledge-network validation, dense correlation workflows,
      prior-guided differential correlation testing, rewiring scores, pathway
      subnetworks, graph creation, and Cytoscape-ready export.

    o Added synthetic example assets, tests, and an introductory vignette.

    o Added permutation validation for node-level rewiring scores.

SIGNIFICANT USER-VISIBLE CHANGES

    o Removed createSparseMultiOmicCorrelations(). Prior-guided sparse
      differential correlation is now performed with
      testDifferentialCorrelation(candidateEdgeTable=...).

    o Adjusted differential-correlation p-values across the full tested
      edge universe before applying significance filters. This behavior was
      changed in version 0.99.1 for the rewiring-first workflow.
