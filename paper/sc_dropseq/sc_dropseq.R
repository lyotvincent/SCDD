library(Seurat)
library(SeuratObject)
library(aricode)
library(ggplot2)
library(ggpubr)
library(ROCR)

# for Umap
label.path <- paste0("data/sc_dropseq.label.txt")
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

sc_dropseq.st <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    st <- getUmap(data, label, has.name, cluster=F)
    return (st)
}

sc_dropseq.cst <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    cst <- getUmap(data, label, has.name, cluster=T)
    return (cst)
}

sc_dropseq.umap <- function(st, title, has.name=F){
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

sc_dropseq.cluster <- function(cst, title, has.name=F){
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
sc_dropseq.cstcomp <- function(cst){
    label <- as.factor(label[,1])
    comp <- clustComp(Idents(cst), label)
    comp <- c(comp, precision.recall(Idents(cst), label))
    return (comp)
}

# paths
raw.path <- paste0("data/sc_dropseq.raw.tsv")
SCDD.path <- paste0("results/SCDD/sc_dropseq_SCDD_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/sc_dropseq_Diffusion_impute.tsv")
scIGANs.path <- paste0("results/scIGANs/sc_dropseq_scIGANs_impute.tsv")
magic.path <- paste0("results/MAGIC/sc_dropseq_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/sc_dropseq_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/sc_dropseq_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/sc_dropseq_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/sc_dropseq_DrImpute_impute.tsv")
scGNN.path <- paste0("results/scGNN/sc_dropseq_scGNN_impute.tsv")
ALRA.path <- paste0("results/ALRA/sc_dropseq_ALRA_impute.tsv")
VIPER.path <- paste0("results/VIPER/sc_dropseq_VIPER_impute.tsv")
scVI.path <- paste0("results/scVI/sc_dropseq_scVI_impute.tsv")

# run Seurat's umap and visualization
raw.st <- sc_dropseq.st(raw.path)
SCDD.st <- sc_dropseq.st(SCDD.path)
Diffusion.st <- sc_dropseq.st(Diffusion.path)
scIGANs.st <- sc_dropseq.st(scIGANs.path)
magic.st <- sc_dropseq.st(magic.path)
saver.st <- sc_dropseq.st(saver.path)
dca.st <- sc_dropseq.st(dca.path)
DeepImpute.st <- sc_dropseq.st(DeepImpute.path)
DrImpute.st <- sc_dropseq.st(DrImpute.path)
scGNN.st <- sc_dropseq.st(scGNN.path)
ALRA.st <- sc_dropseq.st(ALRA.path)
scVI.st <- sc_dropseq.st(scVI.path)
VIPER.st <- sc_dropseq.st(VIPER.path)

raw.umap <- sc_dropseq.umap(raw.st, "Raw")
SCDD.umap <- sc_dropseq.umap(SCDD.st, "SCDD")
Diffusion.umap <- sc_dropseq.umap(Diffusion.st, "SCDD(Diffusion)")
scIGANs.umap <- sc_dropseq.umap(scIGANs.st, "scIGANs")
magic.umap <- sc_dropseq.umap(magic.st, "MAGIC")
saver.umap <- sc_dropseq.umap(saver.st, "SAVER")
dca.umap <- sc_dropseq.umap(dca.st, "DCA")
DeepImpute.umap <- sc_dropseq.umap(DeepImpute.st, "DeepImpute")
DrImpute.umap <- sc_dropseq.umap(DrImpute.st, "DrImpute")
scGNN.umap <- sc_dropseq.umap(scGNN.st, "scGNN")
ALRA.umap <- sc_dropseq.umap(ALRA.st, "ALRA")
scVI.umap <- sc_dropseq.umap(scVI.st, "scVI")
VIPER.umap <- sc_dropseq.umap(VIPER.st, "VIPER")

pdf("paper/sc_dropseq/sc_dropseq_umap31.pdf", 9, 3.5)
ggarrange(raw.umap, SCDD.umap, Diffusion.umap, 
          ncol = 3, nrow = 1, common.legend=T, legend="none")
dev.off()

svg("paper/sc_dropseq/sc_dropseq_umap32.svg", 15, 7.5)
ggarrange(saver.umap, DeepImpute.umap, DrImpute.umap, scIGANs.umap, VIPER.umap,
          scGNN.umap, magic.umap, dca.umap, ALRA.umap, scVI.umap,
          ncol = 5, nrow = 2, common.legend=T, legend="bottom")
dev.off()


# run Seurat's cluster and save the Idents
raw.cst <- sc_dropseq.cst(raw.path)
SCDD.cst <- sc_dropseq.cst(SCDD.path)
Diffusion.cst <- sc_dropseq.cst(Diffusion.path)
magic.cst <- sc_dropseq.cst(magic.path)
saver.cst <- sc_dropseq.cst(saver.path)
dca.cst <- sc_dropseq.cst(dca.path)
DeepImpute.cst <- sc_dropseq.cst(DeepImpute.path)
DrImpute.cst <- sc_dropseq.cst(DrImpute.path)
scGNN.cst <- sc_dropseq.cst(scGNN.path)
ALRA.cst <- sc_dropseq.cst(ALRA.path)
scVI.cst <- sc_dropseq.cst(scVI.path)
VIPER.cst <- sc_dropseq.cst(VIPER.path)

predicts <- data.frame(as.numeric(as.factor(label[,1])),
                 Idents(raw.cst), Idents(SCDD.cst), Idents(Diffusion.cst),
                 Idents(magic.cst), Idents(saver.cst), Idents(dca.cst),
                 Idents(DeepImpute.cst), Idents(DrImpute.cst), Idents(scGNN.cst),
                 Idents(ALRA.cst), Idents(scVI.cst), Idents(VIPER.cst))
write.table(predicts, file='temp/sc_dropseq_predicts1.tsv', sep='\t')

# Visualization of clusters
raw.cluster <- sc_dropseq.cluster(raw.cst, "Raw")
SCDD.cluster <- sc_dropseq.cluster(SCDD.cst, "SCDD")
Diffusion.cluster <- sc_dropseq.cluster(Diffusion.cst, "SCDD(Diffusion)")
magic.cluster <- sc_dropseq.cluster(magic.cst, "MAGIC")
saver.cluster <- sc_dropseq.cluster(saver.cst, "SAVER")
dca.cluster <- sc_dropseq.cluster(dca.cst, "DCA")
DeepImpute.cluster <- sc_dropseq.cluster(DeepImpute.cst, "DeepImpute")
DrImpute.cluster <- sc_dropseq.cluster(DrImpute.cst, "DrImpute")
scGNN.cluster <- sc_dropseq.cluster(scGNN.cst, "scGNN")
ALRA.cluster <- sc_dropseq.cluster(ALRA.cst, "ALRA")
scVI.cluster <- sc_dropseq.cluster(scVI.cst, "scVI")
VIPER.cluster <- sc_dropseq.cluster(VIPER.cst, "VIPER")

pdf("paper/sc_dropseq/sc_dropseq_cluster1.pdf", 12, 12)
ggarrange(raw.cluster, SCDD.cluster, Diffusion.cluster, magic.cluster,
             saver.cluster, dca.cluster, DeepImpute.cluster, DrImpute.cluster,
              scGNN.cluster, ALRA.cluster, scVI.cluster,
              VIPER.cluster, ncol = 3, nrow = 4, common.legend=T, legend="none")
dev.off()

# Competitions among methods
raw.perf <- sc_dropseq.cstcomp(raw.cst)
SCDD.perf <- sc_dropseq.cstcomp(SCDD.cst)
Diffusion.perf <- sc_dropseq.cstcomp(Diffusion.cst)
magic.perf <- sc_dropseq.cstcomp(magic.cst)
saver.perf <- sc_dropseq.cstcomp(saver.cst)
dca.perf <- sc_dropseq.cstcomp(dca.cst)
DeepImpute.perf <- sc_dropseq.cstcomp(DeepImpute.cst)
DrImpute.perf <- sc_dropseq.cstcomp(DrImpute.cst)
scGNN.perf <- sc_dropseq.cstcomp(scGNN.cst)
ALRA.perf <- sc_dropseq.cstcomp(ALRA.cst)
scVI.perf <- sc_dropseq.cstcomp(scVI.cst)
VIPER.perf <- sc_dropseq.cstcomp(VIPER.cst)

perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf,
            magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, scGNN.perf,
            ALRA.perf, scVI.perf, VIPER.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER',
            'DCA', 'DeepImpute', 'DrImpute', 'scGNN', 'ALRA', 'scVI', 'VIPER')
write.table(pf, "temp/sc_dropseq_perfs1.tsv", sep='\t')
