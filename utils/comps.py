import magic
import os
from dca.api import dca
from utils.io import LoadData, SaveData
import rpy2.robjects as robjects
import pandas as pd
import numpy as np
import anndata as ad
import umap
import scanpy as sc
import time

# MAGIC
def imputeByMAGIC(expName, dataPath, labelPath=None):
    modelName = 'MAGIC'
    X, y = LoadData(expName, dataPath, labelPath)
    magic_operator = magic.MAGIC()
    X_magic = magic_operator.fit_transform(X)
    result = X_magic
    result = result.astype(np.int)
    SaveData(expName, modelName, result)
    print("write MAGIC successfully!")

# DCA
def imputeByDCA(expName, dataPath, labelPath=None):
    modelName = 'DCA'
    data = pd.read_csv(dataPath, sep = '\t', index_col = 0)
    _, label = LoadData(expName, dataPath, labelPath)
    data = pd.DataFrame(data.T, dtype = int)
    # obs用于保存细胞的信息
    obs = pd.DataFrame(index=data.index)
    obs['label'] = label
    # vars用于保存基因的信息
    var_names = data.columns
    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(np.array(data), obs=obs, var=var)
    sc.pp.filter_genes(adata, min_counts=1)
    dca(adata, threads=7)
    df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names, dtype=int)
    df = df.T
    if os.path.isdir("results") == False:
        os.makedirs("results")
    dir = "results/"+ modelName
    if os.path.isdir(dir) == False:
        os.makedirs(dir)
    path = dir + "/" + expName + "_" + modelName + "_impute.tsv"
    df.to_csv(path, sep='\t')
    print("write DCA successfully!")

# SAVER
def imputeBySAVER(expName, dataPath, labelPath=None):
    modelName = 'SAVER'
    X, y = LoadData(expName, dataPath, labelPath)
    result = robjects.r('''
    library(SAVER)
    raw <- read.table("%s", header=TRUE, sep="\t")
    rownames(raw) <- raw[,1]
    raw <- raw[, -1]
    saver_res <- saver(raw, ncores=7)
    result <- saver_res$estimate
    return (result)
    ''' % (dataPath))
    result = np.array(result)
    SaveData(expName, modelName, result, needTrans=False)
    print("write SAVER successfully!")