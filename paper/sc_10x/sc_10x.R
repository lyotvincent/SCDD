library(Seurat)
library(SeuratObject)
library(aricode)
library(ggplot2)
library(ggpubr)
library(ROCR)

# for Umap
label.path <- paste0("data/sc_10x.label.txt")
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

sc_10x.st <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    st <- getUmap(data, label, has.name, cluster=F)
    return (st)
}

sc_10x.cst <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    cst <- getUmap(data, label, has.name, cluster=T)
    return (cst)
}

sc_10x.umap <- function(st, title, has.name=F){
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

sc_10x.cluster <- function(cst, title, has.name=F){
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
sc_10x.cstcomp <- function(cst){
    label <- as.factor(label[,1])
    comp <- clustComp(Idents(cst), label)
    comp <- c(comp, precision.recall(Idents(cst), label))
    return (comp)
}

# paths
raw.path <- paste0("data/sc_10x.raw.tsv")
SCDD.path <- paste0("results/SCDD/sc_10x_SCDD1_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/sc_10x_Diffusion_impute.tsv")
scIGANs.path <- paste0("results/scIGANs/sc_10x_scIGANs_impute.tsv")
magic.path <- paste0("results/MAGIC/sc_10x_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/sc_10x_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/sc_10x_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/sc_10x_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/sc_10x_DrImpute_impute.tsv")
scGNN.path <- paste0("results/scGNN/sc_10x_scGNN_impute.tsv")
ALRA.path <- paste0("results/ALRA/sc_10x_ALRA_impute.tsv")
VIPER.path <- paste0("results/VIPER/sc_10x_VIPER_impute.tsv")
scVI.path <- paste0("results/scVI/sc_10x_scVI_impute.tsv")

# run Seurat's umap and visualization
raw.st <- sc_10x.st(raw.path)
SCDD.st <- sc_10x.st(SCDD.path)
Diffusion.st <- sc_10x.st(Diffusion.path)
scIGANs.st <- sc_10x.st(scIGANs.path)
magic.st <- sc_10x.st(magic.path)
saver.st <- sc_10x.st(saver.path)
dca.st <- sc_10x.st(dca.path)
DeepImpute.st <- sc_10x.st(DeepImpute.path)
DrImpute.st <- sc_10x.st(DrImpute.path)
scGNN.st <- sc_10x.st(scGNN.path)
ALRA.st <- sc_10x.st(ALRA.path)
scVI.st <- sc_10x.st(scVI.path)
VIPER.st <- sc_10x.st(VIPER.path)

raw.umap <- sc_10x.umap(raw.st, "Raw")
SCDD.umap <- sc_10x.umap(SCDD.st, "SCDD")
Diffusion.umap <- sc_10x.umap(Diffusion.st, "SCDD(Diffusion)")
scIGANs.umap <- sc_10x.umap(scIGANs.st, "scIGANs")
magic.umap <- sc_10x.umap(magic.st, "MAGIC")
saver.umap <- sc_10x.umap(saver.st, "SAVER")
dca.umap <- sc_10x.umap(dca.st, "DCA")
DeepImpute.umap <- sc_10x.umap(DeepImpute.st, "DeepImpute")
DrImpute.umap <- sc_10x.umap(DrImpute.st, "DrImpute")
scGNN.umap <- sc_10x.umap(scGNN.st, "scGNN")
ALRA.umap <- sc_10x.umap(ALRA.st, "ALRA")
scVI.umap <- sc_10x.umap(scVI.st, "scVI")
VIPER.umap <- sc_10x.umap(VIPER.st, "VIPER")

pdf("paper/sc_10x/sc_10x_umap31.pdf", 9, 3.5)
ggarrange(raw.umap, SCDD.umap, Diffusion.umap, 
          ncol = 3, nrow = 1, common.legend=T, legend="none")
dev.off()

svg("paper/sc_10x/sc_10x_umap32.svg", 15, 7.5)
ggarrange(saver.umap, DeepImpute.umap, DrImpute.umap, scIGANs.umap, VIPER.umap,
          scGNN.umap, magic.umap, dca.umap, ALRA.umap, scVI.umap,
          ncol = 5, nrow = 2, common.legend=T, legend="bottom")
dev.off()


# run Seurat's cluster and save the Idents
raw.cst <- sc_10x.cst(raw.path)
SCDD.cst <- sc_10x.cst(SCDD.path)
Diffusion.cst <- sc_10x.cst(Diffusion.path)
magic.cst <- sc_10x.cst(magic.path)
saver.cst <- sc_10x.cst(saver.path)
dca.cst <- sc_10x.cst(dca.path)
DeepImpute.cst <- sc_10x.cst(DeepImpute.path)
DrImpute.cst <- sc_10x.cst(DrImpute.path)
scGNN.cst <- sc_10x.cst(scGNN.path)
ALRA.cst <- sc_10x.cst(ALRA.path)
scVI.cst <- sc_10x.cst(scVI.path)
VIPER.cst <- sc_10x.cst(VIPER.path)

predicts <- data.frame(as.numeric(as.factor(label[,1])),
                 Idents(raw.cst), Idents(SCDD.cst), Idents(Diffusion.cst),
                 Idents(magic.cst), Idents(saver.cst), Idents(dca.cst),
                 Idents(DeepImpute.cst), Idents(DrImpute.cst), Idents(scGNN.cst),
                 Idents(ALRA.cst), Idents(scVI.cst), Idents(VIPER.cst))
write.table(predicts, file='temp/sc_10x_predicts1.tsv', sep='\t')

# Visualization of clusters
raw.cluster <- sc_10x.cluster(raw.cst, "Raw")
SCDD.cluster <- sc_10x.cluster(SCDD.cst, "SCDD")
Diffusion.cluster <- sc_10x.cluster(Diffusion.cst, "SCDD(Diffusion)")
magic.cluster <- sc_10x.cluster(magic.cst, "MAGIC")
saver.cluster <- sc_10x.cluster(saver.cst, "SAVER")
dca.cluster <- sc_10x.cluster(dca.cst, "DCA")
DeepImpute.cluster <- sc_10x.cluster(DeepImpute.cst, "DeepImpute")
DrImpute.cluster <- sc_10x.cluster(DrImpute.cst, "DrImpute")
scGNN.cluster <- sc_10x.cluster(scGNN.cst, "scGNN")
ALRA.cluster <- sc_10x.cluster(ALRA.cst, "ALRA")
scVI.cluster <- sc_10x.cluster(scVI.cst, "scVI")
VIPER.cluster <- sc_10x.cluster(VIPER.cst, "VIPER")

pdf("paper/sc_10x/sc_10x_cluster1.pdf", 12, 12)
ggarrange(raw.cluster, SCDD.cluster, Diffusion.cluster, magic.cluster,
             saver.cluster, dca.cluster, DeepImpute.cluster, DrImpute.cluster,
              scGNN.cluster, ALRA.cluster, scVI.cluster,
              VIPER.cluster, ncol = 3, nrow = 4, common.legend=T, legend="none")
dev.off()

# Competitions among methods
raw.perf <- sc_10x.cstcomp(raw.cst)
SCDD.perf <- sc_10x.cstcomp(SCDD.cst)
Diffusion.perf <- sc_10x.cstcomp(Diffusion.cst)
magic.perf <- sc_10x.cstcomp(magic.cst)
saver.perf <- sc_10x.cstcomp(saver.cst)
dca.perf <- sc_10x.cstcomp(dca.cst)
DeepImpute.perf <- sc_10x.cstcomp(DeepImpute.cst)
DrImpute.perf <- sc_10x.cstcomp(DrImpute.cst)
scGNN.perf <- sc_10x.cstcomp(scGNN.cst)
ALRA.perf <- sc_10x.cstcomp(ALRA.cst)
scVI.perf <- sc_10x.cstcomp(scVI.cst)
VIPER.perf <- sc_10x.cstcomp(VIPER.cst)

perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf,
            magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, scGNN.perf,
            ALRA.perf, scVI.perf, VIPER.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER',
            'DCA', 'DeepImpute', 'DrImpute', 'scGNN', 'ALRA', 'scVI', 'VIPER')
write.table(pf, "temp/sc_10x_perfs1.tsv", sep='\t')
