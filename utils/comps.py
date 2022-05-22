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
from deepimpute.multinet import MultiNet

# MAGIC
def imputeByMAGIC(expName, dataPath, labelPath=None):
    try:
        modelName = 'MAGIC'
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        magic_operator = magic.MAGIC()
        X_magic = magic_operator.fit_transform(X)
        result = X_magic
        result = result.astype(np.int)
        SaveData(expName, modelName, result)
        print("write MAGIC successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write MAGIC failed!', file=f)

# DCA
def imputeByDCA(expName, dataPath, labelPath=None):
    try:
        modelName = 'DCA'
        data = pd.read_csv(dataPath, sep = '\t', index_col = 0)
        _, label = LoadData(expName, dataPath, labelPath=labelPath)
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
    except:
        with open('log.txt', 'a+') as f:
            print('write DCA failed!', file=f)

# SAVER
def imputeBySAVER(expName, dataPath, labelPath=None):
    try:
        modelName = 'SAVER'
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
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
    except:
        with open('log.txt', 'a+') as f:
            print('write SAVER failed!', file=f)

# DrImpute
def imputeByDrImpute(expName, dataPath, labelPath=None):
    modelName = 'DrImpute'
    X, y = LoadData(expName, dataPath, labelPath=labelPath)
    try:
        result = robjects.r('''
            library(DrImpute)
            library(SummarizedExperiment)
            library(scDatasets)
            library(Matrix)
            raw <- read.table("%s", header=TRUE, sep="\t")
            rownames(raw) <- raw[,1]
            raw <- raw[, -1]
            raw <- Matrix(as.matrix(raw))
            raw.log <- log(raw + 1)
            set.seed(1)
            result <- DrImpute(raw.log, mc.cores = 7)
            result <- exp(result) - 1
            return(result)
            ''' % (dataPath))
        result = np.array(result)
        SaveData(expName, modelName, result, needTrans=False)
        print("write DrImpute successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write DrImpute failed!', file=f)

# VIPER
def imputeByVIPER(expName, dataPath, labelPath=None):
    try:
        modelName = "VIPER"
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        result = robjects.r('''
            library(VIPER)
            raw <- read.table("%s", header=TRUE, sep="\t")
            rownames(raw) <- raw[,1]
            raw <- raw[, -1]
            res <- VIPER(raw, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                         report = FALSE, outdir = NULL, prefix = NULL)
            return (res$imputed)
            ''' % (dataPath))
        result = np.array(result)
        SaveData(expName, modelName, result, needTrans=False)
        print("write VIPER successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write VIPER failed!', file=f)

# scIGANs
def imputeByscIGANs(expName, dataPath, labelPath=None):
    try:
        modelName = "scIGANs"
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        cmd = "scIGANs/scIGANs {0} -b {1} -o results/scIGANs -j {2}".format(dataPath, 128, expName)
        os.system(cmd)
        print("write scIGANs successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write scIGANs failed!', file=f)

# DeepImpute
def imputeByDeepImpute(expName, dataPath, labelPath=None):
    try:
        modelName = "DeepImpute"
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        data = pd.DataFrame(X)
        model = MultiNet()
        model.fit(data)
        result = model.predict(data)
        SaveData(expName, modelName, result)
        print("write DeepImpute successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write DeepImpute failed!', file=f)

# scTSSR
def imputeByscTSSR(expName, dataPath, labelPath=None):
    modelName = "scTSSR"
    X, y = LoadData(expName, dataPath, labelPath=labelPath)
    result = robjects.r('''
        library(scTSSR)
        raw <- read.table("%s", header=TRUE, sep="\t")
        rownames(raw) <- raw[,1]
        raw <- raw[, -1]
        result <- scTSSR(raw, 
                          percent=0.05, 
                          learning_rate=0.0001, 
                          epochs=100, 
                          run_batch = FALSE,
                          estimates.only = TRUE,
                          ncores = 1)
        return (result)
        ''' % (dataPath))
    result = np.array(result)
    SaveData(expName, modelName, result, needTrans=False)
    print("write scTSSR successfully!")

# scImpute
def imputeByscImpute(expName, dataPath, labelPath=None):
    try:
        modelName = "scImpute"
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        result = robjects.r('''
            library(scImpute)
            res <- scimpute(# full path to raw count matrix
                  count_path = "%s", 
                  infile = "txt",           # format of input file
                  outfile = "txt",          # format of output file
                  out_dir = "./results/scImpute/",           # full path to output directory
                  labeled = FALSE,          # cell type labels not available
                  Kcluster = 6,
                  drop_thre = 0.5,          # threshold set on dropout probability
                  ncores = 1)              # number of cores used in parallel computation
            ''' % (dataPath))
        tpfiles = os.listdir("results/scImpute")
        tpfiles.remove("scimpute_count.txt")
        for file in tpfiles:
            os.remove("results/scImpute/" + file)
        os.rename("results/scImpute/scimpute_count.txt",
                  "results/scImpute/{0}_{1}_impute.tsv".format(expName, modelName))
        print("write scImpute successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write scImpute failed!', file=f)

# scGNN
def imputeByscGNN(expName, dataPath, labelPath=None):
    try:
        modelName = "scGNN"
        if not os.path.isdir("csvdata"):
            os.makedirs("csvdata")
        datasetName = expName + ".raw.csv"
        csvdataPath = "csvdata/" + datasetName
        df = pd.read_csv(dataPath, sep='\t', index_col=0)
        df.to_csv(csvdataPath, header=True)
        pp_cmd = "~/anaconda3/bin/python " \
                 "-W ignore scGNN/PreprocessingscGNN.py " \
                 "--datasetName {0} " \
                 "--datasetDir csvdata/ " \
                 "--LTMGDir csvdata/ " \
                 "--filetype CSV " \
                 "--geneSelectnum 15000 " \
                 "--inferLTMGTag".format(datasetName)
        os.system(pp_cmd)
        ipt_cmd = " ~/anaconda3/bin/python " \
                  "-W ignore scGNN/scGNN.py " \
                  "--datasetName csvdata " \
                  "--datasetDir ./ " \
                  "--LTMGDir ./ " \
                  "--outputDir results/scGNN/ " \
                  "--EM-iteration 2 " \
                  "--Regu-epochs 50 " \
                  "--EM-epochs 20 " \
                  "--quickmode " \
                  "--nonsparseMode " \
                  "--regulized-type LTMG"
        os.system(ipt_cmd)
        df = pd.read_csv("results/scGNN/csvdata_recon.csv", index_col=0)
        df = np.exp(df) - 1
        tpfile = os.listdir("results/scGNN")
        for f in tpfile:
            os.remove("results/scGNN/" + f)
        df.to_csv("results/scGNN/{0}_{1}_impute.tsv".format(expName, modelName), sep='\t', header=True)
        csvfile = os.listdir("csvdata")
        for f in csvfile:
            os.remove("csvdata/" + f)
        os.rmdir("csvdata")
        print("write scGNN successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write scGNN failed!', file=f)

# EnImpute
def imputeByEnImpute(expName, dataPath, labelPath=None):
    modelName = "scImpute"
    X, y = LoadData(expName, dataPath, labelPath=labelPath)


# ALRA
def imputeByALRA(expName, dataPath, labelPath=None):
    modelName = "ALRA"
    try:
        X, y = LoadData(expName, dataPath, labelPath=labelPath)
        result = robjects.r('''
                source('ALRA/alra.R')
                raw <- read.table("%s", header=TRUE, sep='\t')
                rownames(raw) <- raw[,1]
                raw <- raw[, -1]
                count <- t(as.matrix(raw))
                A_norm <- normalize_data(count)
                k_choice <- choose_k(A_norm)
                A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]
                return(A_norm_completed)
                    ''' % (dataPath))
        result = np.array(result)
        SaveData(expName, modelName, result, needTrans=True)
        print("write ALRA successfully!")
    except:
        with open('log.txt', 'a+') as f:
            print('write ALRA failed!', file=f)





