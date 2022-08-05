library(Seurat)
library(SeuratObject)
library(aricode)
library(ROCR)

getUMAP <- function(raw, label, has_name = TRUE) {
    if (has_name == FALSE) {
        rownames(raw) <- raw[, 1] # 将第一列作为行名
        raw <- raw[, -1] # 去除第一列
    }
    dt <- which(label == "Unknown") # 找到标签为Unknown的数据后删除
    label <- label[-dt, ]
    raw <- raw[, -dt]
    raw <- data.matrix(raw) # 转化成matrix形式
    dimnames <- list(rownames(raw), colnames(raw)) # 保存行名和列名
    Seurat_raw <- CreateSeuratObject(raw, project = "SeuratObj", assay = "RNA")
    Seurat_raw@active.ident <- as.factor(label)
    Seurat_raw <- NormalizeData(Seurat_raw, normalization.method = "LogNormalize", scale.factor = 10000)
    Seurat_raw <- FindVariableFeatures(Seurat_raw, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(Seurat_raw)
    Seurat_raw <- ScaleData(Seurat_raw, features = all.genes)
    Seurat_raw <- RunUMAP(Seurat_raw, features = VariableFeatures(object = Seurat_raw))
    return(Seurat_raw)
}

cluster <- function(raw, label, has_name = TRUE, do_norm = TRUE) {
    if (has_name == FALSE) {
        rownames(raw) <- raw[, 1] # 将第一列作为行名
        raw <- raw[, -1] # 去除第一列
    }
    raw <- data.matrix(raw) # 转化成matrix形式
    dimnames <- list(rownames(raw), colnames(raw)) # 保存行名和列名
    Seurat_raw <- CreateSeuratObject(raw, project = "SeuratObj", assay = "RNA")
    raw <- data.matrix(raw) # 转化成matrix形式
    dimnames <- list(rownames(raw), colnames(raw)) # 保存行名和列名
    Seurat_raw <- CreateSeuratObject(raw, project = "SeuratObj", assay = "RNA")
    Seurat_raw@meta.data$orig.ident <- as.factor(label[, 1])
    if(do_norm){
      Seurat_raw <- NormalizeData(Seurat_raw, normalization.method = "LogNormalize", scale.factor = 10000)
    }
    Seurat_raw <- FindVariableFeatures(Seurat_raw, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(Seurat_raw)
    Seurat_raw <- ScaleData(Seurat_raw, features = all.genes)
    Seurat_raw <- RunPCA(Seurat_raw, features = VariableFeatures(object = Seurat_raw))
    Seurat_raw <- FindNeighbors(Seurat_raw, dims = 1:6)
    Seurat_raw <- FindClusters(Seurat_raw, resolution = 0.2)
    Seurat_raw <- RunUMAP(Seurat_raw, features = VariableFeatures(object = Seurat_raw))
    return(Seurat_raw)
}

no_uk_ident <- function(st, label) {
    dt <- which(label == "Unknown")
    ident <- Idents(st)
    ident <- ident[-dt]
    return(ident)
}

no_uk_comp <- function(st, label) {
    ident <- no_uk_ident(st, label)
    dt <- which(label == "Unknown")
    ulabel <- label[-dt, ]
    Comp <- clustComp(ident, as.factor(ulabel))
    perf <- c(Comp, precision.recall(ident, as.factor(ulabel)))
    return(perf)
}

## defined a function to compute the precision, recall, and F value of two cluster results.
precision.recall <- function(c1, c2) {
    ## c1, a vector of reference labels that assigned to each cluster
    ## c2, a vector of predicted labels that assigned to each cluster
    # c1=true.label
    # c2=km.label.k
    ref.label <- as.vector(outer(c1, c1, "=="))
    pred.label <- as.vector(outer(c2, c2, "=="))
    pred <- prediction(as.numeric(pred.label), ref.label)

    # Accuracy
    acc.tmp <- performance(pred, "acc")
    acc <- as.numeric(acc.tmp@y.values[[1]][2])

    # ROC curve
    ROC.perf <- performance(pred, "tpr", "fpr")
    # ROC area under the curve
    auc.tmp <- performance(pred, "auc")
    auc <- as.numeric(auc.tmp@y.values)

    ## precision
    prec.tmp <- performance(pred, "ppv")
    prec <- as.numeric(prec.tmp@y.values[[1]][2])
    ## F1-score
    f.tmp <- performance(pred, "f")
    f <- as.numeric(f.tmp@y.values[[1]][2])
    ##
    return(list(F.score = f, AUC = auc, ACC = acc))
}

raw.path <- paste0("data/Bone.raw.tsv")
label.path <- paste0("data/Bone.label.txt")
SCDD.path <- paste0("results/SCDD/Bone_SCDD_impute.tsv")
Diffusion.path <- paste0("results/Diffusion/Bone_Diffusion_impute.tsv")
magic.path <- paste0("results/MAGIC/Bone_MAGIC_impute.tsv")
saver.path <- paste0("results/SAVER/Bone_SAVER_impute.tsv")
dca.path <- paste0("results/DCA/Bone_DCA_impute.tsv")
DeepImpute.path <- paste0("results/DeepImpute/Bone_DeepImpute_impute.tsv")
DrImpute.path <- paste0("results/DrImpute/Bone_DrImpute_impute.tsv")
scGNN.path <- paste0("results/scGNN/Bone_scGNN_impute.tsv")
ALRA.path <- paste0("results/ALRA/Bone_ALRA_impute.tsv")
scIGANs.path <- paste0("results/scIGANs/Bone_scIGANs_impute.tsv")
scVI.path <- paste0("results/scVI/Bone_scVI_impute.tsv")
scTSSR.path <- paste0("results/scTSSR/Bone_scTSSR_impute.tsv")
EnImpute.path <- paste0("results/EnImpute/Bone_EnImpute_impute.tsv")

raw <- read.table(raw.path, sep = '\t',  header = TRUE)
label <- read.table(label.path, sep = '\t',  header = FALSE)
SCDD <- read.table(SCDD.path, sep = '\t',  header = TRUE)
Diffusion <- read.table(Diffusion.path, sep = '\t',  header = TRUE)
magic <- read.table(magic.path, sep = '\t',  header = TRUE)
saver <- read.table(saver.path, sep = '\t',  header = TRUE)
dca <- read.table(dca.path, sep = '\t',  header = TRUE)
DeepImpute <- read.table(DeepImpute.path, sep = '\t',  header = TRUE)
DrImpute <- read.table(DrImpute.path, sep = '\t',  header = TRUE)
scGNN <- read.table(scGNN.path, sep = '\t',  header = TRUE)
ALRA <- read.table(ALRA.path, sep = '\t',  header = TRUE)
scIGANs <- read.table(scIGANs.path, sep = '\t',  header = TRUE)
scVI <- read.table(scVI.path, sep = '\t',  header = TRUE)
scTSSR <- read.table(scTSSR.path, sep = '\t',  header = TRUE)
EnImpute <- read.table(EnImpute.path, sep = '\t',  header = TRUE)


raw.cluster <- cluster(raw, label, F)
SCDD.cluster <- cluster(SCDD, label, F)
Diffusion.cluster <- cluster(Diffusion, label, F)
magic.cluster <- cluster(magic, label, F)
saver.cluster <- cluster(saver, label, T)
dca.cluster <- cluster(dca, label, F)
DeepImpute.cluster <- cluster(DeepImpute, label, F)
DrImpute.cluster <- cluster(DrImpute, label, F)
scGNN.cluster <- cluster(scGNN, label, F)
ALRA.cluster <- cluster(ALRA, label, F, do_norm = FALSE)
scIGANs.cluster <- cluster(scIGANs, label, T)
scVI.cluster <- cluster(scVI, label, F)
scTSSR.cluster <- cluster(scTSSR, label, T)
EnImpute.cluster <- cluster(EnImpute, label, T)

raw.perf <- no_uk_comp(raw.cluster, label)
SCDD.perf <- no_uk_comp(SCDD.cluster, label)
Diffusion.perf <- no_uk_comp(Diffusion.cluster, label)
magic.perf <- no_uk_comp(magic.cluster, label)
saver.perf <- no_uk_comp(saver.cluster, label)
dca.perf <- no_uk_comp(dca.cluster, label)
DeepImpute.perf <- no_uk_comp(DeepImpute.cluster, label)
DrImpute.perf <- no_uk_comp(DrImpute.cluster, label)
scGNN.perf <- no_uk_comp(scGNN.cluster, label)
ALRA.perf <- no_uk_comp(ALRA.cluster, label)
scIGANs.perf <- no_uk_comp(scIGANs.cluster, label)
scVI.perf <- no_uk_comp(scVI.cluster, label)
scTSSR.perf <- no_uk_comp(scTSSR.cluster, label)
EnImpute.perf <- no_uk_comp(EnImpute.cluster, label)


perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf, magic.perf, saver.perf, dca.perf,
            DeepImpute.perf, DrImpute.perf, scGNN.perf, ALRA.perf, scIGANs.perf, scVI.perf, scTSSR.perf, EnImpute.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER', 'DCA',
                  'DeepImpute', 'DrImpute', 'scGNN', 'ALRA', 'scIGANs', 'scVI', 'scTSSR', 'EnImpute')
# perfs1 <- read.table("temp/Bone_perfs1.tsv")
write.table(pf, "temp/Bone_perfs1.tsv", sep='\t')
# perfs1 <- rbind(perfs1, scTSSR.perf, EnImpute.perf)

raw.ident <- no_uk_ident(raw.cluster, label)
SCDD.ident <- no_uk_ident(SCDD.cluster, label)
Diffusion.ident <- no_uk_ident(Diffusion.cluster, label)
magic.ident <- no_uk_ident(magic.cluster, label)
saver.ident <- no_uk_ident(saver.cluster, label)
dca.ident <- no_uk_ident(dca.cluster, label)
DeepImpute.ident <- no_uk_ident(DeepImpute.cluster, label)
DrImpute.ident <- no_uk_ident(DrImpute.cluster, label)
scGNN.ident <- no_uk_ident(scGNN.cluster, label)
ALRA.ident <- no_uk_ident(ALRA.cluster, label)
scIGANs.ident <- no_uk_ident(scIGANs.cluster, label)
scVI.ident <- no_uk_ident(scVI.cluster, label)
scTSSR.ident <- no_uk_ident(scTSSR.cluster, label)
EnImpute.ident <- no_uk_ident(EnImpute.cluster, label)

dt <- which(label=='Unknown')    # 找到标签为Unknown的数据后删除
ulabel <- label[-dt,]
predicts <- data.frame(as.numeric(as.factor(ulabel)), # save the idents to cal purity
                  raw.ident, SCDD.ident, Diffusion.ident,
                  magic.ident, saver.ident, dca.ident,
                  DeepImpute.ident, DrImpute.ident, scGNN.ident,
                  ALRA.ident, scIGANs.ident, scVI.ident)
write.table(predicts, file='temp/Bone_predicts1.tsv', sep='\t')


# raw.cluster <- cluster(raw, label, F)
# raw.perf <- no_uk_comp(raw.cluster, label)
# perf <- data.frame(raw.perf)
# rownames(perf) <- c('raw')

# raw.ident <- no_uk_ident(raw.cluster, label)
# dt <- which(label == 'Unknown')    # 找到标签为Unknown的数据后删除
# ulabel <- label[-dt, ]
# res <- data.frame(as.numeric(as.factor(ulabel)), # save the idents to cal purity
#                   raw.ident)

# write.table(res, file='temp/Bone_predicts.tsv', sep='\t')
