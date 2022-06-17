import numpy as np
from numpy import random
import pandas as pd
import random
from tqdm import tqdm
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix, save_npz, load_npz
import os

cells = np.array([])
genes = np.array([])


def LoadData(expName, dataPath, format="tsv",
             labelPath=None, needTrans=True, labelHeader=None,
             structure='Array', slot="raw"):
    """
    :rtype: tuple(np.array, np.array)
    :param expName: experiment name
    :param dataPath: the origin data path
    :param labelPath: the label path
    :param needTrans: if true, transpose the origin data.
    :param labelHeader: need to read the first line or not.
    :param structure: data sturcture used in our method, Array or AnnData.
    :return: training data and training label
    """
    global cells
    global genes
    print("Experiment:{0}".format(expName))
    print("Data path:{0}".format(dataPath))
    print("Label path:{0}".format(labelPath))
    print("Using Structure:{0}".format(structure))

    # Load labels
    label = list()
    if labelPath:
        train_label = np.array(pd.read_csv(labelPath, header=labelHeader))
        train_label = [elem[0] for elem in train_label]
        tp, l = dict(), 0
        cell_type = list()
        for elem in train_label:
            if elem in tp.keys():
                label.append(tp[elem])
            else:
                cell_type.append(elem)
                tp[elem] = l
                label.append(l)
                l += 1

    # Load data for 2 sturctures form tsv or h5ad
    if structure == "AnnData":
        train_data = ad.AnnData()
        if format == "h5ad":
            adata = ad.read_h5ad(dataPath)
            cells = adata.obs_names
            genes = adata.var_names
            if slot == 'raw':
                train_data = ad.AnnData(X=adata.raw.X,
                              obs=pd.DataFrame(index=adata.obs_names),
                              var=pd.DataFrame(index=adata.var_names))
            else:
                train_data = adata
            train_data.obs['total_counts'] = train_data.X.sum(axis=1).A1
            train_data.var['total_counts'] = train_data.X.sum(axis=0).A1
        elif format == "tsv":
            adata = sc.read_csv(dataPath, delimiter='\t', dtype='int')
            if needTrans:
                adata = adata.transpose()
            cells = adata.obs_names
            genes = adata.var_names
            train_data = adata
            print("Data size:{0}".format(adata.shape))
            train_data.obs['total_counts'] = train_data.X.sum(axis=1)
            train_data.var['total_counts'] = train_data.X.sum(axis=0)
        print("Load Data to AnnData OK.")
        return train_data, label
    else:
        df = pd.DataFrame()
        if format == "tsv":
            df = pd.read_csv(dataPath, sep='\t', index_col=0)
            cells = df.columns
            genes = df.index
        elif format == "h5ad":
            adata = ad.read_h5ad(dataPath)
            df = pd.DataFrame(adata.raw.X.todense(),
                          index=adata.obs_names,
                          columns=adata.var_names,
                          dtype='int')
            cells = df.index
            genes = df.columns
        # set global values to save cells and genes, which will be used in SaveData
        if needTrans and format != 'h5ad':
            train_data = np.array(df).transpose()
        else:
            train_data = np.array(df)
        print("Data size:{0}".format(train_data.shape))
        print("Load Data to Array OK.")
        return train_data, label


def SaveData(expName, modelName, res, format="tsv", needTrans=True, id = None):
    """
    This function will call `impute' member in model and save the result,
        which will use global cells and genes.
    Please make sure that genes are in index and cells are in rows.
    We will save the result in `results' dictionary
    :param needTrans:  if true, transpose the origin data.
    :param expName:  the experiment name
    :param modelName:  the model name
    :param res:  impute result
    """
    if os.path.isdir("results") == False:
        os.makedirs("results")
    dir = "results/"+ modelName
    if os.path.isdir(dir) == False:
        os.makedirs(dir)
    ID = str("")
    if id:
        ID = str(id)
    if format == "tsv":
        path = dir + "/" + expName + "_" + modelName + ID + "_impute.tsv"
        if needTrans:
            res = pd.DataFrame(res.transpose())
        else:
            res = pd.DataFrame(res)
        res.index = genes
        res.columns = cells
        res.to_csv(path, sep='\t')
        print("Write to {0} successfully!".format(path))
    elif format == "h5ad":
        path = dir + "/" + expName + "_" + modelName + ID + "_impute.h5ad"
        adata = ad.AnnData(X=res, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
        adata.write_h5ad(path)
        print("Write to {0} successfully!".format(path))


def SaveTargets(M, Omega, Target, dropout_rate, null_genes, sparse=False):
    if sparse == False:
        np.savez("temp/clusters_M", M)
    else:
        save_npz("temp/clusters_M.npz", M)
    np.savez("temp/Omega", Omega)
    np.savez("temp/dropout_rate", dropout_rate)
    np.savez("temp/null_genes", null_genes)
    np.savez("temp/Target", Target)

def LoadTargets(sparse=False):
    if sparse == False:
        M = np.load("temp/clusters_M.npz")
        M = M['arr_0']
    else:
        M = load_npz("temp/clusters_M.npz")
    print("Load related Matrix `M' from temp...")
    Omega = np.load("temp/Omega.npz")
    Omega = Omega['arr_0']
    print("Load Omega from temp...")
    Target = np.load("temp/Target.npz")
    Target = Target['arr_0']
    print("Load Target from temp...")
    dropout_rate = np.load("temp/dropout_rate.npz")
    dropout_rate = dropout_rate['arr_0']
    print("Load dropout_rate from temp...")
    null_genes = np.load("temp/null_genes.npz")
    null_genes = null_genes['arr_0']
    print("Load null_genes from temp...")
    return M, Omega, Target, dropout_rate, null_genes

def DataSplit(h5adpath, outdir, chunksize=5000):
    if os.path.isdir(outdir) == False:
        os.makedirs(outdir)
    import anndata as ad
    adata = ad.read_h5ad(h5adpath)
    adsize = adata.shape[0]
    adidx = list(range(adsize))
    random.shuffle(adidx)
    partnum = int((adsize-1) / chunksize + 1)
    for i in tqdm(range(partnum)):
        end = chunksize * (i+1)
        if end > adsize:
            end = adsize
        idx = adidx[chunksize * i : end]
        adata0 = adata[idx, :]
        adata0.obs['idx'] = idx
        adata0.write_h5ad(outdir + '/chunk{0}_part{1}.h5ad'.format(chunksize, i))
    print("data splited to {0}".format("outdir"))
    # df = pd.DataFrame(adata0.T.X.todense(),
    #                   index=adata0.obs_names,
    #                   columns=adata0.var_names,
    #                   dtype='int')
    # lbl = adata0.obs['free_annotation']
    # df.to_csv(outdir + '/cnt/chunk{0}_part{1}.tsv'.format(chunksize, i), sep = '\t', header=True, chunksize=1000)
    # lbl.to_csv(outdir + '/lbl/chunk{0}_part{1}.txt'.format(chunksize, i), index=False, header=False)
    # idx = pd.DataFrame(idx)
    # idx.to_csv(outdir + '/idx/chunk{0}_part{1}.txt'.format(chunksize, i), index=False, header=False)

def DataMerge(datadir, imputedir):
    datas = os.listdir(datadir)

    imputed = os.listdir(imputedir)
    part = ad.read_h5ad(imputedir + imputed[0])
    data = ad.read_h5ad(datadir + datas[0])
    part.obs = data.obs
    for i in range(1, len(imputed)):
        print("datai:{0};imputedi:{1}".format(datas[i], imputed[i]))
        parti = ad.read_h5ad(imputedir + imputed[i])
        datai = ad.read_h5ad(datadir + datas[i])
        parti.obs = datai.obs
        part = ad.concat([part, parti])
    part.write_h5ad(imputedir + "immune_impute_all.h5ad")

def DataMask(dataPath, mask=0.1):
    df = pd.read_csv(dataPath, sep='\t', index_col=0)
    rand = np.random.uniform(-mask, 1 - mask, df.shape) > 0
    mask_df = np.where(rand, df, np.zeros(df.shape))
    mask_df = pd.DataFrame(mask_df, index=df.index, columns=df.columns)
    mask_pt = pd.DataFrame(rand, index=df.index, columns=df.columns)
    mask_df.to_csv("data/mask{0}_".format(mask) + dataPath.split('/')[1], sep='\t', header=True)
    mask_pt.to_csv("data/maskpt{0}_".format(mask) + dataPath.split('/')[1], sep='\t', header=True)

def DataMask2(dataPath, p=0.8):
    df = pd.read_csv(dataPath, sep='\t', index_col=0)
    df = df.T
    gene_depth = np.sum(df, axis=0)
    cell_depth = np.sum(df, axis=1)
    gene_median = np.median(gene_depth)
    cell_median = np.median(cell_depth)
    gene = gene_depth > gene_median
    cell = cell_depth > cell_median
    df0 = df[cell]
    df0 = df0.T
    df1 = df0[gene]
    A0 = df1.T
    A1 = A0.apply(lambda x:random.binomial(n=x, p=p, size=144))
    A0.to_csv("data/DS_Cellcycle.raw.txt", sep='\t', header=True)
    A1.to_csv("data/DS0.8_Cellcycle.raw.txt", sep='\t', header=True)





