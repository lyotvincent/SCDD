library(Seurat)
library(SeuratObject)
library(aricode)
library(ggplot2)
library(ggpubr)
library(ROCR)

# for Umap
label.path <- paste0("data/sc_celseq2.label.txt")
label <- read.table(label.path, sep = '\t',  header = FALSE)


getUmap <- function(data, label, has.name = F, cluster=F){
  if(has.name == FALSE){
    rownames(data) <- data[, 1]    # 将第一列作为行名
    data <- data[, -1]    # 去除第一列
  }
  data <- data.matrix(data)    # 转化成matrix形式
  dimnames = list(rownames(data), colnames(data))    # 保存行名和列名
  st <- CreateSeuratObject(data, project = "Seurat", assay = "RNA")
  st@meta.data$orig.ident = label[,1]
  st@active.ident = as.factor(label[,1])
  st <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000)
  st <- FindVariableFeatures(st, selection.method = "vst", nfeatures = 4000)
  all.genes <- rownames(st)
  st <- ScaleData(st, features = all.genes)
  st <- RunUMAP(st, features = VariableFeatures(object = st))
  if(cluster==T){
      st <- RunPCA(st, features = VariableFeatures(object = st))
      st<- FindNeighbors(st, dims = 1:6)
      st <- FindClusters(st, resolution = 0.2)
  }
  return(st)
}

sc_celseq2.st <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    st <- getUmap(data, label, has.name, cluster=F)
    return (st)
}

sc_celseq2.cst <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    cst <- getUmap(data, label, has.name, cluster=T)
    return (cst)
}

sc_celseq2.umap <- function(st, title, has.name=F){
    umap.plot <- DimPlot(st, reduction = "umap", pt.size = 0.8)+
      labs(x="UMAP 1",y="UMAP 2",title=title)+
      theme(title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face="plain"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
    return (umap.plot)
}

sc_celseq2.cluster <- function(cst, title, has.name=F){
    umap.plot <- DimPlot(cst, reduction = "umap", pt.size = 0.8)+
      labs(x="UMAP 1",y="UMAP 2",title=title)+
      theme(title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face="plain"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
    return (umap.plot)
}

## defined a function to compute the precision, recall, and F value of two cluster results.
precision.recall <- function(c1, c2){
  ##c1, a vector of reference labels that assigned to each cluster
  ##c2, a vector of predicted labels that assigned to each cluster
  #c1=true.label
  #c2=km.label.k
  ref.label <- as.vector(outer(c1, c1, "=="))
  pred.label <- as.vector(outer(c2, c2, "=="))
  pred <- prediction(as.numeric(pred.label),ref.label)

  # Accuracy
  acc.tmp <- performance(pred, "acc");
  acc <- as.numeric(acc.tmp@y.values[[1]][2])

  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr");
  # ROC area under the curve
  auc.tmp <- performance(pred,"auc");
  auc <- as.numeric(auc.tmp@y.values)

  ##precision
  prec.tmp = performance(pred, "ppv")
  prec <- as.numeric(prec.tmp@y.values[[1]][2])
  ## F1-score
  f.tmp = performance(pred, "f")
  f <- as.numeric(f.tmp@y.values[[1]][2])
  return(list(F.score=f, AUC=auc, ACC=acc))
}

# Performances among methods
sc_celseq2.cstcomp <- function(cst){
    label <- as.factor(label[,1])
    comp <- clustComp(Idents(cst), label)
    comp <- c(comp, precision.recall(Idents(cst), label))
    return (comp)
}

# paths
raw.path <- paste0("data/sc_celseq2.raw.tsv")
SCDD.path <- paste0("results/SCDD/sc_celseq2_SCDD1_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/sc_celseq2_Diffusion_impute.tsv")
magic.path <- paste0("results/MAGIC/sc_celseq2_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/sc_celseq2_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/sc_celseq2_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/sc_celseq2_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/sc_celseq2_DrImpute_impute.tsv")
scGNN.path <- paste0("results/scGNN/sc_celseq2_scGNN_impute.tsv")
ALRA.path <- paste0("results/ALRA/sc_celseq2_ALRA_impute.tsv")
VIPER.path <- paste0("results/VIPER/sc_celseq2_VIPER_impute.tsv")
scVI.path <- paste0("results/scVI/sc_celseq2_scVI_impute.tsv")

# run Seurat's umap and visualization
raw.st <- sc_celseq2.st(raw.path)
SCDD.st <- sc_celseq2.st(SCDD.path)
Diffusion.st <- sc_celseq2.st(Diffusion.path)
magic.st <- sc_celseq2.st(magic.path)
saver.st <- sc_celseq2.st(saver.path)
dca.st <- sc_celseq2.st(dca.path)
DeepImpute.st <- sc_celseq2.st(DeepImpute.path)
DrImpute.st <- sc_celseq2.st(DrImpute.path)
scGNN.st <- sc_celseq2.st(scGNN.path)
ALRA.st <- sc_celseq2.st(ALRA.path)
scVI.st <- sc_celseq2.st(scVI.path)
VIPER.st <- sc_celseq2.st(VIPER.path)

raw.umap <- sc_celseq2.umap(raw.st, "Raw")
SCDD.umap <- sc_celseq2.umap(SCDD.st, "SCDD")
Diffusion.umap <- sc_celseq2.umap(Diffusion.st, "SCDD(Diffusion)")
magic.umap <- sc_celseq2.umap(magic.st, "MAGIC")
saver.umap <- sc_celseq2.umap(saver.st, "SAVER")
dca.umap <- sc_celseq2.umap(dca.st, "DCA")
DeepImpute.umap <- sc_celseq2.umap(DeepImpute.st, "DeepImpute")
DrImpute.umap <- sc_celseq2.umap(DrImpute.st, "DrImpute")
scGNN.umap <- sc_celseq2.umap(scGNN.st, "scGNN")
ALRA.umap <- sc_celseq2.umap(ALRA.st, "ALRA")
scVI.umap <- sc_celseq2.umap(scVI.st, "scVI")
VIPER.umap <- sc_celseq2.umap(VIPER.st, "VIPER")

pdf("paper/sc_celseq2/sc_celseq2_umap1.pdf", 12, 12)
ggarrange(raw.umap, SCDD.umap, Diffusion.umap, magic.umap,
             saver.umap, dca.umap, DeepImpute.umap, DrImpute.umap,
              scGNN.umap, ALRA.umap, scVI.umap, VIPER.umap,
              ncol = 3, nrow = 4, common.legend=T, legend="bottom")
dev.off()

# run Seurat's cluster and save the Idents
raw.cst <- sc_celseq2.cst(raw.path)
SCDD.cst <- sc_celseq2.cst(SCDD.path)
Diffusion.cst <- sc_celseq2.cst(Diffusion.path)
magic.cst <- sc_celseq2.cst(magic.path)
saver.cst <- sc_celseq2.cst(saver.path)
dca.cst <- sc_celseq2.cst(dca.path)
DeepImpute.cst <- sc_celseq2.cst(DeepImpute.path)
DrImpute.cst <- sc_celseq2.cst(DrImpute.path)
scGNN.cst <- sc_celseq2.cst(scGNN.path)
ALRA.cst <- sc_celseq2.cst(ALRA.path)
scVI.cst <- sc_celseq2.cst(scVI.path)
VIPER.cst <- sc_celseq2.cst(VIPER.path)

predicts <- data.frame(as.numeric(as.factor(label[,1])),
                 Idents(raw.cst), Idents(SCDD.cst), Idents(Diffusion.cst),
                 Idents(magic.cst), Idents(saver.cst), Idents(dca.cst),
                 Idents(DeepImpute.cst), Idents(DrImpute.cst), Idents(scGNN.cst),
                 Idents(ALRA.cst), Idents(scVI.cst), Idents(VIPER.cst))
write.table(predicts, file='temp/sc_celseq2_predicts1.tsv', sep='\t')

# Visualization of clusters
raw.cluster <- sc_celseq2.cluster(raw.cst, "Raw")
SCDD.cluster <- sc_celseq2.cluster(SCDD.cst, "SCDD")
Diffusion.cluster <- sc_celseq2.cluster(Diffusion.cst, "SCDD(Diffusion)")
magic.cluster <- sc_celseq2.cluster(magic.cst, "MAGIC")
saver.cluster <- sc_celseq2.cluster(saver.cst, "SAVER")
dca.cluster <- sc_celseq2.cluster(dca.cst, "DCA")
DeepImpute.cluster <- sc_celseq2.cluster(DeepImpute.cst, "DeepImpute")
DrImpute.cluster <- sc_celseq2.cluster(DrImpute.cst, "DrImpute")
scGNN.cluster <- sc_celseq2.cluster(scGNN.cst, "scGNN")
ALRA.cluster <- sc_celseq2.cluster(ALRA.cst, "ALRA")
scVI.cluster <- sc_celseq2.cluster(scVI.cst, "scVI")
VIPER.cluster <- sc_celseq2.cluster(VIPER.cst, "VIPER")

pdf("paper/sc_celseq2/sc_celseq2_cluster1.pdf", 12, 12)
ggarrange(raw.cluster, SCDD.cluster, Diffusion.cluster, magic.cluster,
             saver.cluster, dca.cluster, DeepImpute.cluster, DrImpute.cluster,
              scGNN.cluster, ALRA.cluster, scVI.cluster,
              VIPER.cluster, ncol = 3, nrow = 4, common.legend=T, legend="none")
dev.off()

# Competitions among methods
raw.perf <- sc_celseq2.cstcomp(raw.cst)
SCDD.perf <- sc_celseq2.cstcomp(SCDD.cst)
Diffusion.perf <- sc_celseq2.cstcomp(Diffusion.cst)
magic.perf <- sc_celseq2.cstcomp(magic.cst)
saver.perf <- sc_celseq2.cstcomp(saver.cst)
dca.perf <- sc_celseq2.cstcomp(dca.cst)
DeepImpute.perf <- sc_celseq2.cstcomp(DeepImpute.cst)
DrImpute.perf <- sc_celseq2.cstcomp(DrImpute.cst)
scGNN.perf <- sc_celseq2.cstcomp(scGNN.cst)
ALRA.perf <- sc_celseq2.cstcomp(ALRA.cst)
scVI.perf <- sc_celseq2.cstcomp(scVI.cst)
VIPER.perf <- sc_celseq2.cstcomp(VIPER.cst)

perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf,
            magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, scGNN.perf,
            ALRA.perf, scVI.perf, VIPER.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER',
            'DCA', 'DeepImpute', 'DrImpute', 'scGNN', 'ALRA', 'scVI', 'VIPER')
write.table(pf, "temp/sc_celseq2_perfs1.tsv", sep='\t')
