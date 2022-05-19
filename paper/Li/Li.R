library(Seurat)
library(SeuratObject)
library(aricode)
library(ggplot2)
library(ggpubr)
library(ROCR)

# for Umap
label.path <- paste0("data/Li.label.txt")
label <- read.table(label.path, sep = '\t',  header = FALSE)
# for Comps
cluster.path <- paste0("data/Li.cluster.txt")
cluster <- read.table(cluster.path, sep = '\t',  header = FALSE)


getUmap <- function(data, label, has.name = F, cluster=F){
  if(has.name == FALSE){
    rownames(data) <- data[, 1]    # 将第一列作为行名
    data <- data[, -1]    # 去除第一列
  }
  data <- data.matrix(data)    # 转化成matrix形式
  dimnames = list(rownames(data), colnames(data))    # 保存行名和列名
  st <- CreateSeuratObject(data, project = "StomachSeuratObj", assay = "RNA")
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

Li.st <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    st <- getUmap(data, label, has.name, cluster=F)
    return (st)
}

Li.cst <- function(path, has.name=F){
    data <- read.table(path, sep = '\t',  header = TRUE)
    cst <- getUmap(data, label, has.name, cluster=T)
    return (cst)
}

Li.umap <- function(st, title, has.name=F){
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

Li.cluster <- function(cst, title, has.name=F){
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
Li.cstcomp <- function(cst){
    label <- as.factor(cluster[,1])
    comp <- clustComp(Idents(cst), label)
    comp <- c(comp, precision.recall(Idents(cst), label))
    return (comp)
}

# paths
raw.path <- paste0("data/Li.raw.tsv")
SCDD.path <- paste0("results/SCDD/Li_SCDD1_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/Li_Diffusion_impute.tsv")
magic.path <- paste0("results/MAGIC/Li_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/Li_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/Li_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/Li_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/Li_DrImpute_impute.tsv")
VIPER.path <- paste0("results/VIPER/Li_VIPER_impute.tsv")
scGNN.path <- paste0("results/DeepImpute/Li_scGNN_impute.tsv")
scImpute.path <- paste0("results/DrImpute/Li_scImpute_impute.tsv")
scTSSR.path <- paste0("results/scTSSR/Li_scTSSR_impute.tsv")

# run Seurat's umap and visualization
raw.st <- Li.st(raw.path)
SCDD.st <- Li.st(SCDD.path)
Diffusion.st <- Li.st(Diffusion.path)
magic.st <- Li.st(magic.path)
saver.st <- Li.st(saver.path)
dca.st <- Li.st(dca.path)
DeepImpute.st <- Li.st(DeepImpute.path)
DrImpute.st <- Li.st(DrImpute.path)
VIPER.st <- Li.st(VIPER.path)
scGNN.st <- Li.st(scGNN.path)
scImpute.st <- Li.st(scImpute.path)
scTSSR.st <- Li.st(scTSSR.path)

raw.umap <- Li.umap(raw.st, "Raw")
SCDD.umap <- Li.umap(SCDD.st, "SCDD")
Diffusion.umap <- Li.umap(Diffusion.st, "SCDD(Diffusion)")
magic.umap <- Li.umap(magic.st, "MAGIC")
saver.umap <- Li.umap(saver.st, "SAVER")
dca.umap <- Li.umap(dca.st, "DCA")
DeepImpute.umap <- Li.umap(DeepImpute.st, "DeepImpute")
DrImpute.umap <- Li.umap(DrImpute.st, "DrImpute")
VIPER.umap <- Li.umap(VIPER.st, "VIPER")
scGNN.umap <- Li.umap(scGNN.st, "scGNN")
scImpute.umap <- Li.umap(scImpute.st, "scImpute")
scTSSR.umap <- Li.umap(scTSSR.st, "scTSSR")

pdf("paper/Li/Li_umap1.pdf", 12, 12)
ggarrange(raw.umap, SCDD.umap, Diffusion.umap, magic.umap, 
             saver.umap, dca.umap, DeepImpute.umap, DrImpute.umap,
              VIPER.umap, scGNN.umap, scImpute.umap, scTSSR.umap,
              ncol = 3, nrow = 2, common.legend=T, legend="bottom")
dev.off()

# run Seurat's cluster and save the Idents
raw.cst <- Li.cst(raw.path)
SCDD.cst <- Li.cst(SCDD.path)
Diffusion.cst <- Li.cst(Diffusion.path)
magic.cst <- Li.cst(magic.path)
saver.cst <- Li.cst(saver.path)
dca.cst <- Li.cst(dca.path)
DeepImpute.cst <- Li.cst(DeepImpute.path)
DrImpute.cst <- Li.cst(DrImpute.path)
VIPER.cst <- Li.cst(VIPER.path)
scGNN.cst <- Li.cst(scGNN.path)
scImpute.cst <- Li.cst(scImpute.path)
scTSSR.cst <- Li.cst(scTSSR.path)

predicts <- data.frame(as.numeric(as.factor(cluster[,1])), 
                 Idents(raw.cst), Idents(SCDD.cst), Idents(Diffusion.cst), 
                 Idents(magic.cst), Idents(saver.cst), Idents(dca.cst),
                 Idents(DeepImpute.cst), Idents(DrImpute.cst), Idents(VIPER.cst),
                 Idents(scGNN.cst), Idents(scImpute.cst), Idents(scTSSR.cst))
write.table(predicts, file='temp/Li_predicts1.tsv', sep='\t')

# Visualization of clusters
raw.cluster <- Li.cluster(raw.cst, "Raw")
SCDD.cluster <- Li.cluster(SCDD.cst, "SCDD")
Diffusion.cluster <- Li.cluster(Diffusion.cst, "SCDD(Diffusion)")
magic.cluster <- Li.cluster(magic.cst, "MAGIC")
saver.cluster <- Li.cluster(saver.cst, "SAVER")
dca.cluster <- Li.cluster(dca.cst, "DCA")
DeepImpute.cluster <- Li.cluster(DeepImpute.cst, "DeepImpute")
DrImpute.cluster <- Li.cluster(DrImpute.cst, "DrImpute")
VIPER.cluster <- Li.cluster(VIPER.cst, "VIPER")
scGNN.cluster <- Li.cluster(scGNN.cst, "scGNN")
scImpute.cluster <- Li.cluster(scImpute.cst, "scImpute")
scTSSR.cluster <- Li.cluster(scTSSR.cst, "scTSSR")

pdf("paper/Li/Li_cluster1.pdf", 12, 12)
ggarrange(raw.cluster, SCDD.cluster, Diffusion.cluster, magic.cluster, 
             saver.cluster, dca.cluster, DeepImpute.cluster, DrImpute.cluster,
              VIPER.cluster, scGNN.cluster, scImpute.cluster,
              scTSSR.cluster, ncol = 3, nrow = 2, common.legend=T, legend="none")
dev.off()

# Competitions among methods
raw.perf <- Li.cstcomp(raw.cst)
SCDD.perf <- Li.cstcomp(SCDD.cst)
Diffusion.perf <- Li.cstcomp(Diffusion.cst)
magic.perf <- Li.cstcomp(magic.cst)
saver.perf <- Li.cstcomp(saver.cst)
dca.perf <- Li.cstcomp(dca.cst)
DeepImpute.perf <- Li.cstcomp(DeepImpute.cst)
DrImpute.perf <- Li.cstcomp(DrImpute.cst)
VIPER.perf <- Li.cstcomp(VIPER.cst)
scGNN.perf <- Li.cstcomp(scGNN.cst)
scImpute.perf <- Li.cstcomp(scImpute.cst)
scTSSR.perf <- Li.cstcomp(scTSSR.cst)

perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf,
            magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, VIPER.perf,
            scGNN.perf, scImpute.perf, scTSSR.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER',
            'DCA', 'DeepImpute', 'DrImpute', 'VIPER', 'scGNN', 'scImpute', 'scTSSR')
write.table(pf, "temp/Li_perfs1.tsv", sep='\t')
