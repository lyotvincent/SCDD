library(Seurat)
library(SeuratObject)
library(aricode)
library(ggplot2)
library(ggpubr)
library(ROCR)

# for Umap
label.path <- paste0("data/CIDR.label.txt")
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
      st <- FindClusters(st, resolution = 0.5)
  }
  return(st)
}

CIDR.st <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    st <- getUmap(data, label, has.name, cluster=F)
    return (st)
}

CIDR.cst <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    cst <- getUmap(data, label, has.name, cluster=T)
    return (cst)
}

CIDR.umap <- function(st, title, has.name=F){
    umap.plot <- DimPlot(st, reduction = "umap", pt.size = 0.8)+
      labs(x="UMAP 1",y="UMAP 2",title=title)+
      theme(title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face="plain"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
    return (umap.plot)
}

CIDR.cluster <- function(cst, title, has.name=F){
    umap.plot <- DimPlot(cst, reduction = "umap", pt.size = 0.8)+
      labs(x="UMAP 1",y="UMAP 2",title=title)+
      theme(title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face="plain"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines"),
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
CIDR.cstcomp <- function(cst){
    label <- as.factor(label[,1])
    comp <- clustComp(Idents(cst), label)
    comp <- c(comp, precision.recall(Idents(cst), label))
    return (comp)
}

# paths
raw.path <- paste0("data/CIDR.raw.tsv")
SCDD.path <- paste0("results/SCDD/CIDR_SCDD1_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/CIDR_Diffusion_impute.tsv")
scIGANs.path <- paste0("results/scIGANs/CIDR_scIGANs_impute.tsv")
magic.path <- paste0("results/MAGIC/CIDR_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/CIDR_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/CIDR_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/CIDR_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/CIDR_DrImpute_impute.tsv")
scGNN.path <- paste0("results/scGNN/CIDR_scGNN_impute.tsv")
ALRA.path <- paste0("results/ALRA/CIDR_ALRA_impute.tsv")
VIPER.path <- paste0("results/VIPER/CIDR_VIPER_impute.tsv")
scVI.path <- paste0("results/scVI/CIDR_scVI_impute.tsv")

# run Seurat's umap and visualization
raw.st <- CIDR.st(raw.path)
SCDD.st <- CIDR.st(SCDD.path)
scIGANs.st <- CIDR.st(scIGANs.path)
Diffusion.st <- CIDR.st(Diffusion.path)
magic.st <- CIDR.st(magic.path)
saver.st <- CIDR.st(saver.path)
dca.st <- CIDR.st(dca.path)
DeepImpute.st <- CIDR.st(DeepImpute.path)
DrImpute.st <- CIDR.st(DrImpute.path)
scGNN.st <- CIDR.st(scGNN.path)
ALRA.st <- CIDR.st(ALRA.path)
scVI.st <- CIDR.st(scVI.path)
VIPER.st <- CIDR.st(VIPER.path)

raw.umap <- CIDR.umap(raw.st, "Raw")
SCDD.umap <- CIDR.umap(SCDD.st, "SCDD")
Diffusion.umap <- CIDR.umap(Diffusion.st, "SCDD(Diffusion)")
scIGANs.umap <- CIDR.umap(Diffusion.st, "scIGANs")
magic.umap <- CIDR.umap(magic.st, "MAGIC")
saver.umap <- CIDR.umap(saver.st, "SAVER")
dca.umap <- CIDR.umap(dca.st, "DCA")
DeepImpute.umap <- CIDR.umap(DeepImpute.st, "DeepImpute")
DrImpute.umap <- CIDR.umap(DrImpute.st, "DrImpute")
scGNN.umap <- CIDR.umap(scGNN.st, "scGNN")
ALRA.umap <- CIDR.umap(ALRA.st, "ALRA")
scVI.umap <- CIDR.umap(scVI.st, "scVI")
VIPER.umap <- CIDR.umap(VIPER.st, "VIPER")

svg("paper/CIDR/CIDR_umap31.svg", 6.5, 2.4)
ggarrange(raw.umap, SCDD.umap, Diffusion.umap, 
          ncol = 3, nrow = 1, common.legend=T, legend="none")
dev.off()

svg("paper/CIDR/CIDR_umap32.svg", 10, 4.6)
ggarrange(saver.umap, DeepImpute.umap, DrImpute.umap, scIGANs.umap, VIPER.umap,
          scGNN.umap, magic.umap, dca.umap, ALRA.umap, scVI.umap,
          ncol = 5, nrow = 2, common.legend=T, legend="bottom")
dev.off()

# run Seurat's cluster and save the Idents
raw.cst <- CIDR.cst(raw.path)
SCDD.cst <- CIDR.cst(SCDD.path)
Diffusion.cst <- CIDR.cst(Diffusion.path)
magic.cst <- CIDR.cst(magic.path)
saver.cst <- CIDR.cst(saver.path)
dca.cst <- CIDR.cst(dca.path)
DeepImpute.cst <- CIDR.cst(DeepImpute.path)
DrImpute.cst <- CIDR.cst(DrImpute.path)
scGNN.cst <- CIDR.cst(scGNN.path)
ALRA.cst <- CIDR.cst(ALRA.path)
scVI.cst <- CIDR.cst(scVI.path)
VIPER.cst <- CIDR.cst(VIPER.path)

predicts <- data.frame(as.numeric(as.factor(label[,1])),
                 Idents(raw.cst), Idents(SCDD.cst), Idents(Diffusion.cst),
                 Idents(magic.cst), Idents(saver.cst), Idents(dca.cst),
                 Idents(DeepImpute.cst), Idents(DrImpute.cst), Idents(scGNN.cst),
                 Idents(ALRA.cst), Idents(scVI.cst), Idents(VIPER.cst))
write.table(predicts, file='temp/CIDR_predicts1.tsv', sep='\t')

# Visualization of clusters
raw.cluster <- CIDR.cluster(raw.cst, "Raw")
SCDD.cluster <- CIDR.cluster(SCDD.cst, "SCDD")
Diffusion.cluster <- CIDR.cluster(Diffusion.cst, "SCDD(Diffusion)")
scIGANs.cluster <- CIDR.cluster(Diffusion.cst, "SCDD(Diffusion)")
magic.cluster <- CIDR.cluster(magic.cst, "MAGIC")
saver.cluster <- CIDR.cluster(saver.cst, "SAVER")
dca.cluster <- CIDR.cluster(dca.cst, "DCA")
DeepImpute.cluster <- CIDR.cluster(DeepImpute.cst, "DeepImpute")
DrImpute.cluster <- CIDR.cluster(DrImpute.cst, "DrImpute")
scGNN.cluster <- CIDR.cluster(scGNN.cst, "scGNN")
ALRA.cluster <- CIDR.cluster(ALRA.cst, "ALRA")
scVI.cluster <- CIDR.cluster(scVI.cst, "scVI")
VIPER.cluster <- CIDR.cluster(VIPER.cst, "VIPER")

svg("paper/CIDR/CIDR_cluster1.svg", 12, 12)
ggarrange(raw.cluster, SCDD.cluster, Diffusion.cluster, magic.cluster,
             saver.cluster, dca.cluster, DeepImpute.cluster, DrImpute.cluster,
              scGNN.cluster, ALRA.cluster, scVI.cluster,
              VIPER.cluster, ncol = 3, nrow = 4, common.legend=T, legend="none")
dev.off()

# Competitions among methods
raw.perf <- CIDR.cstcomp(raw.cst)
SCDD.perf <- CIDR.cstcomp(SCDD.cst)
Diffusion.perf <- CIDR.cstcomp(Diffusion.cst)
magic.perf <- CIDR.cstcomp(magic.cst)
saver.perf <- CIDR.cstcomp(saver.cst)
dca.perf <- CIDR.cstcomp(dca.cst)
DeepImpute.perf <- CIDR.cstcomp(DeepImpute.cst)
DrImpute.perf <- CIDR.cstcomp(DrImpute.cst)
scGNN.perf <- CIDR.cstcomp(scGNN.cst)
ALRA.perf <- CIDR.cstcomp(ALRA.cst)
scVI.perf <- CIDR.cstcomp(scVI.cst)
VIPER.perf <- CIDR.cstcomp(VIPER.cst)

perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf,
            magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, scGNN.perf,
            ALRA.perf, scVI.perf, VIPER.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER',
            'DCA', 'DeepImpute', 'DrImpute', 'scGNN', 'ALRA', 'scVI', 'VIPER')
write.table(pf, "temp/CIDR_perfs1.tsv", sep='\t')
