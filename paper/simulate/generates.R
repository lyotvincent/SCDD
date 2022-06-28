# using splatter to generate scalable dataset
library('splatter')
params <- newSplatParams()
params.groups <- newSplatParams(batchCells = 1000000, nGenes = 1000)
sim <- splatSimulateGroups(params.groups, group.prob = c(0.3, 0.3,0.4),de.prob = 0.4,verbose = TRUE)
write.table(assays(sim)$TrueCounts,file = "data/test_truedata.tsv", sep='\t', quote=F)
write.table(assays(sim)$counts,file = "data/test_data.tsv", sep='\t', quote=F)
write.table(colData(sim)$Group,file = "data/test_label.txt", sep='\t', quote=F)