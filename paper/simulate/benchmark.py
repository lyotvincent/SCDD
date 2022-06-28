# please run this file under the main dictionary
import anndata as ad
import scanpy as sc
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt

dataPath_100w_20clu = ['data/test_data.tsv',
                       'results/SCDD/test_impute.tsv',
                       'results/DCA/test_DCA_impute.tsv',
                       'results/DeepImpute/test_DeepImpute_impute.tsv']

# we test the deep-learing based method in big dataset to test the scalability of SCDD
mtName_100w_20clu = ['Raw', 'SCDD', 'DCA', 'DeepImpute']

label = "data/test_label.txt"

# using ari to benchmark the scalable dataset
def benchmark_ari(dataPath, mtName, labelPath):
    adata = ad.read_csv(dataPath, delimiter='\t')
    adata = adata.T
    # qadata = ad.AnnData(adata.raw.X, obs=adata.obs, var=adata.var)
    label = pd.read_csv(labelPath, sep='\t')
    label.index = adata.obs_names
    adata.obs['free_annotation'] = label['x']
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
    print("Calculating leiden...")
    sc.tl.leiden(adata, key_added='clusters')
    ARI = metrics.adjusted_rand_score(adata.obs['clusters'], adata.obs['free_annotation'])
    # adata.obs['clusters'].to_csv("temp/test100w_SCDD10_clusters.tsv", sep='\t', index=True, header=True)
    print("{0} imputed:{1}".format(mtName, ARI))


for (dt, m) in zip(dataPath_100w_20clu, mtName_100w_20clu):
    benchmark_ari(dt, m, label)
