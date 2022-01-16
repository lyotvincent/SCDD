library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)

cc.pp <- function(st_data, st_label, st_genes, all=T, has.name=F, Ercc=F){
    colnames(st_genes) <- c("ID", "gene")
    st_genes[which(st_genes$gene == '-'), 'gene'] <- st_genes[which(st_genes$gene == '-'), 'ID']
    st_genes[which(duplicated(st_genes$gene)), 'gene'] <- st_genes[which(duplicated(st_genes$gene)), 'ID']
    rownames(st_genes) <- st_genes$ID
    if(all == F){
        rownames(st_data) <- lapply(st_data[,1], function(x){x <- st_genes[x, 'gene']})
    }
    else {
        rownames(st_data) <- st_genes[, 'gene']
    }
    if(has.name == F){
        st_data <- st_data[, -1]    # 去除第一列 
    }
    if(Ercc==F){
        for(elem in c("Ercc1", "Ercc6l2", "Ercc8", "Ercc2", "Ercc3", "Ercc4", "Ercc5", "Ercc6", "Ercc6l")){
            st_data <- st_data[which(rownames(st_data) != elem),]
        }
    }
    return (st_data)
}

cc.plot <- function(st_data, st_label, title, Ercc=F){
    st <- data.matrix(st_data)    # 转化成matrix形式
    dimnames = list(rownames(st), colnames(st))    # 保存行名和列名
    Seurat_st <- CreateSeuratObject(st, project = "SeuratObj", assay = "RNA")
    Seurat_st@meta.data$orig.ident = st_label[,1]
    Idents(Seurat_st) = st_label[,1]
    Seurat_st <- NormalizeData(Seurat_st, normalization.method = "LogNormalize", scale.factor = 10000)
    all.genes <- rownames(Seurat_st)
    Seurat_st <- ScaleData(Seurat_st, features = all.genes)
    if(Ercc == F){
        Seurat_st <- RunUMAP(Seurat_st, features = all.genes)
    }
    else{
        Seurat_st <- RunUMAP(Seurat_st, 
            features = as.factor(c("Ercc1", "Ercc6l2", "Ercc8", 
                                   "Ercc2", "Ercc3", "Ercc4", 
                                   "Ercc5", "Ercc6", "Ercc6l")))
    }
    cluster.plot <- DimPlot(Seurat_st, reduction = "umap", pt.size = 0.8)+ 
      labs(x="UMAP 1",y="UMAP 2",title=title) + 
      theme(plot.title = element_text(hjust = 0.5, face="plain", size = 12),
            legend.title=element_text("Celltypes"),
            legend.position="bottom",
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 12),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA))
    return (cluster.plot)
}


noErcc <- function(path, title, all=T, has.name=F) {
    data <- read.table(path, sep='\t', header=T)
    ne_data <- cc.pp(data, st_label, st_genes, all=all, has.name=has.name, Ercc=F)
    ne_plot <- cc.plot(ne_data, st_label, title, Ercc=F)
    return (ne_plot)
}

Ercc <- function(path, title, all=T, has.name=F) {
    data <- read.table(path, sep='\t', header=T)
    e_data <- cc.pp(data, st_label, st_genes, all=all, has.name=has.name, Ercc=T)
    e_plot <- cc.plot(e_data, st_label, title, Ercc=T)
    return (e_plot)
}

st_label <- read.table("data/Cellcycle.label.txt", sep = '\t',  header = FALSE)
st_genes <- read.table("data/ID.map.txt", sep = '\t',  header = FALSE)
st_genes <- st_genes[, 1:2]

raw.path <- "data/Cellcycle.raw.txt"
SCDD.path <- "results/SCDD/Cellcycle_SCDD_impute.tsv"
Diffusion.path <- "results/Diffusion/Cellcycle_Diffusion_impute.tsv"
magic.path <- "results/MAGIC/Cellcycle_MAGIC_impute.tsv"
saver.path <- "results/SAVER/Cellcycle_SAVER_impute.tsv"
dca.path <- "results/DCA/Cellcycle_DCA_impute.tsv"

# total plot without Erccs
raw.ne_plot <- noErcc(raw.path, "Raw")
SCDD.ne_plot <- noErcc(SCDD.path, "SCDD")
Diffusion.ne_plot <- noErcc(Diffusion.path, "SCDD(Diffusion)")
magic.ne_plot <- noErcc(magic.path, "MAGIC")
saver.ne_plot <- noErcc(saver.path, "SAVER")
dca.ne_plot <- noErcc(dca.path, "DCA", all=F)

# Ercc plot
raw.e_plot <- Ercc(raw.path, "Raw")
SCDD.e_plot <- Ercc(SCDD.path, "SCDD")
Diffusion.e_plot <- Ercc(Diffusion.path, "SCDD(Diffusion)")
magic.e_plot <- Ercc(magic.path, "MAGIC")
saver.e_plot <- Ercc(saver.path, "SAVER")
dca.e_plot <- Ercc(dca.path, "DCA", all=F)

noErcc.plot.s <- ggarrange(raw.ne_plot, SCDD.ne_plot, Diffusion.ne_plot, 
                       magic.ne_plot, saver.ne_plot, dca.ne_plot,
             nrow=2, ncol=3, legend="bottom", common.legend=T)

Ercc.plot.s <- ggarrange(raw.e_plot, SCDD.e_plot, Diffusion.e_plot,
                       magic.e_plot, saver.e_plot, dca.e_plot,
             nrow=2, ncol=3, legend="bottom", common.legend=T)

# Figures
pdf("paper/Cellcycle/Cellcycle_noErcc_s.pdf", 9, 6)
noErcc.plot.s
dev.off()

pdf("paper/Cellcycle/Cellcycle_Ercc_s.pdf", 9, 6)
Ercc.plot.s
dev.off()