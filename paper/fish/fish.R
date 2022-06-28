# Generating gini index, running this script may cost a long time,
# please checkout the results from paper/fish/fish.gini.tsv.
library(reldist)
library(dplyr)
library(data.table)

fs <- read.table("data/fish.subset.txt")
# calculate the gini index for fish data
fs.gi <- apply(fs, 2, function(x){
  return(gini(na.omit(x)))
})
genes <- colnames(fs)
# fish <- read.table("paper/fish/fish0.gini.tsv")

cal_gini <- function(fish){
  fish <- data.frame(fish)
  rownames(fish) <- fish[, 1]
  fish <- fish[, -1]
  fish <- t(fish[genes, ])
  gi <- apply(fish, 2, gini)
  return(gi)
}

raw.fish <- fread("data/fish.raw.tsv", sep='\t')
raw.gi <- cal_gini(raw.fish)

# calculate the gini index for imputed data
SCDD.fish <- fread("results/SCDD/fish_SCDD7_impute.tsv")
SCDD.gi <- cal_gini(SCDD.fish)

SAVER.fish <- fread("results/SAVER/fish_SAVER_impute.tsv")
SAVER.gi <- cal_gini(SAVER.fish)
rm(SAVER.fish)

DeepImpute.fish <- fread("results/DeepImpute/fish_DeepImpute_impute.tsv")
DeepImpute.gi <- cal_gini(DeepImpute.fish)
rm(DeepImpute.fish)

MAGIC.fish <- fread("results/MAGIC/fish_MAGIC_impute.tsv")
MAGIC.gi <- cal_gini(MAGIC.fish)
rm(MAGIC.fish)

DrImpute.fish <- fread("results/DrImpute/fish_DrImpute_impute.tsv")
DrImpute.gi <- cal_gini(DrImpute.fish)
rm(DrImpute.fish)

scGNN.fish <- fread("results/scGNN/fish_scGNN_impute.tsv")
scGNN.gi <- cal_gini(scGNN.fish)
rm(scGNN.fish)

ALRA.fish <- fread("results/ALRA/fish_ALRA_impute.tsv")
ALRA.gi <- cal_gini(ALRA.fish)
rm(ALRA.fish)

DCA.fish <- fread("results/DCA/fish_DCA_impute.tsv")
DCA.gi <- cal_gini(DCA.fish)
rm(DCA.fish)

fish0 <- cbind(fs.gi, raw.gi, SCDD.gi, SAVER.gi, DeepImpute.gi, MAGIC.gi,
               DrImpute.gi, scGNN.gi, ALRA.gi, DCA.gi)

write.table(fish0, "paper/fish/fish0.gini.tsv", sep='\t', quote=F)
