# please run this file under the main dictionary
# this dataset has 260k+ cells with 60k genes
import anndata as ad
import scanpy as sc
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt

dataPath_100w_20clu = ['data/TS_immune.h5ad',
                       'results/SCDD/immune_SCDD_impute.h5ad',
                       'results/DCA/immune_DCA_impute.h5ad',
                       'results/DeepImpute/immune_DeepImpute_impute.h5ad']

# we test the deep-learing based method in big dataset to test the scalability of SCDD
mtName_100w_20clu = ['Raw', 'SCDD', 'DCA', 'DeepImpute']

raw = ad.read_h5ad('data/TS_immune.h5ad')
label = raw.obs['free_annotation']


# using ari to benchmark the scalable dataset
def benchmark_ari(dataPath, mtName, slot=""):
    adata = ad.read_h5ad(dataPath)
    if slot == 'raw':
        adata = ad.AnnData(adata.raw.X, obs=adata.obs, var=adata.var)
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=5)
    sc.tl.leiden(adata, key_added='clusters')
    ARI = metrics.adjusted_rand_score(adata.obs['clusters'], label)
    print("{0} ARI:{1}".format(mtName, ARI))
    # adata.obs['clusters'].to_csv("temp/immune_raw_clusters.tsv", sep='\t', index=True, header=True)


for (dt, m) in zip(dataPath_100w_20clu, mtName_100w_20clu):
    if m == 'Raw':
        benchmark_ari(dt, m, 'raw')
    else:
        benchmark_ari(dt, m)
