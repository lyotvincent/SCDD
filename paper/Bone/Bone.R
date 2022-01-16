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

cluster <- function(raw, label, has_name = TRUE) {
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
    Seurat_raw <- NormalizeData(Seurat_raw, normalization.method = "LogNormalize", scale.factor = 10000)
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

raw <- read.table(raw.path, sep = '\t',  header = TRUE)
label <- read.table(label.path, sep = '\t',  header = FALSE)
SCDD <- read.table(SCDD.path, sep = '\t',  header = TRUE)
Diffusion <- read.table(Diffusion.path, sep = '\t',  header = TRUE)
magic <- read.table(magic.path, sep = '\t',  header = TRUE)
saver <- read.table(saver.path, sep = '\t',  header = TRUE)
dca <- read.table(dca.path, sep = '\t',  header = TRUE)


raw.cluster <- cluster(raw, label, F)
SCDD.cluster <- cluster(SCDD, label, F)
Diffusion.cluster <- cluster(Diffusion, label, F)
magic.cluster <- cluster(magic, label, F)
saver.cluster <- cluster(saver, label, T)
dca.cluster <- cluster(dca, label, F)

raw.perf <- no_uk_comp(raw.cluster, label)
SCDD.perf <- no_uk_comp(SCDD.cluster, label)
Diffusion.perf <- no_uk_comp(Diffusion.cluster, label)
magic.perf <- no_uk_comp(magic.cluster, label)
saver.perf <- no_uk_comp(saver.cluster, label)
dca.perf <- no_uk_comp(dca.cluster, label)
perf <- data.frame(raw.perf)
pf <- rbind(perf, SCDD.perf, Diffusion.perf, magic.perf, saver.perf, dca.perf)
rownames(pf) <- c('Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC', 'SAVER', 'DCA')
write.table(pf, "temp/Bone_perfs.tsv", sep='\t')


raw.ident <- no_uk_ident(raw.cluster, label)
SCDD.ident <- no_uk_ident(SCDD.cluster, label)
Diffusion.ident <- no_uk_ident(Diffusion.cluster, label)
magic.ident <- no_uk_ident(magic.cluster, label)
saver.ident <- no_uk_ident(saver.cluster, label)
dca.ident <- no_uk_ident(dca.cluster, label)
dt <- which(label=='Unknown')    # 找到标签为Unknown的数据后删除
ulabel <- label[-dt,]
predicts <- data.frame(as.numeric(as.factor(ulabel)), # save the idents to cal purity
                  raw.ident, SCDD.ident, Diffusion.ident,
                  magic.ident, saver.ident, dca.ident)
write.table(predicts, file='temp/Bone_predicts.tsv', sep='\t')


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