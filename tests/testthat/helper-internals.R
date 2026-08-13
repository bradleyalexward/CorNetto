## Internal helpers exercised directly by the unit tests. testthat sources
## helper-*.R before the test files, so these are available everywhere and each
## internal is resolved once rather than at every use site.

coerceEdgeTable <- getFromNamespace(".coerceStandardEdgeTable", "CorNetto")
computeDegreeMatched <- getFromNamespace(".computeDegreeMatchedScores", "CorNetto")
computePairwise <- getFromNamespace(".computeCorrelationByPairwiseTests", "CorNetto")
computePearson <- getFromNamespace(".computePearsonCorrelations", "CorNetto")
duplicateUndirected <- getFromNamespace(".duplicateUndirectedEdges", "CorNetto")
finalizeResults <- getFromNamespace(".finalizeDifferentialCorrelationResults", "CorNetto")
fisherDifference <- getFromNamespace(".fisherCorrelationDifference", "CorNetto")
getNodeTable <- getFromNamespace(".getStoredNodeTable", "CorNetto")
makeJobs <- getFromNamespace(".makePermutationJobs", "CorNetto")
nodeKey <- getFromNamespace(".nodeKey", "CorNetto")
permuteLabels <- getFromNamespace(".permuteGroupLabels", "CorNetto")
resolveRewiringPlotLabels <- getFromNamespace(".resolveRewiringPlotLabels", "CorNetto")
setNodeTable <- getFromNamespace(".setStoredNodeTable", "CorNetto")
standardColumns <- getFromNamespace(".standardEdgeColumns", "CorNetto")()
tailProbability <- getFromNamespace(".permutationTailProbability", "CorNetto")
