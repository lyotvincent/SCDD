from SCDD_api import *

if __name__ == '__main__':
    neighbors = 20
    e = SCDD(name="Li", raw="data/Li.raw.tsv", id=980, format="tsv", batch_size=700, structure="AnnData", neighbor_method="SNN")
    e.run(store=False)

    # import anndata as ad
    # import scanpy as sc
    # from sklearn import metrics
    # adata = ad.read_h5ad("data/TS_immune.h5ad")
    # adata = ad.AnnData(adata.raw.X, obs=adata.obs, var=adata.var)
    # sc.pp.filter_genes(adata, min_counts=1)
    # sc.pp.normalize_per_cell(adata)
    # sc.pp.log1p(adata)
    # sc.pp.scale(adata)
    # sc.tl.pca(adata, svd_solver='arpack')
    # sc.pp.neighbors(adata, n_neighbors=50, n_pcs=5)
    # sc.tl.leiden(adata, key_added='clusters')
    # ARI =  metrics.adjusted_rand_score(adata.obs['clusters'], adata.obs['free_annotation'])
    # print("raw:{0}".format(ARI))
    # adata0 = ad.read_h5ad("results/SCDD/immune_SCDD2_impute.h5ad")
    # sc.pp.filter_genes(adata0, min_counts=1)
    # sc.pp.normalize_per_cell(adata0)
    # sc.pp.log1p(adata0)
    # sc.pp.scale(adata0)
    # sc.tl.pca(adata0, svd_solver='arpack')
    # sc.pp.neighbors(adata0, n_neighbors=50, n_pcs=5)
    # sc.tl.leiden(adata0, key_added='clusters')
    # ARI =  metrics.adjusted_rand_score(adata0.obs['clusters'], adata.obs['free_annotation'])
    # print("imputed:{0}".format(ARI))