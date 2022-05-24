import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import anndata as ad
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

s_dataPath = ['data/liver.raw.tsv', 
            'results/SCDD/liver_SCDD_impute.tsv',
            'results/Diffusion/liver_Diffusion_impute.tsv',
            'results/MAGIC/liver_MAGIC_impute.tsv',
            'results/SAVER/liver_SAVER_impute.tsv',
            'results/DCA/liver_DCA_impute.tsv',
            'results/scIGANs/liver_scIGANs_impute.tsv',
            'results/DrImpute/liver_DrImpute_impute.tsv',
            'results/DeepImpute/liver_DeepImpute_impute.tsv',
            'results/scGNN/liver_scGNN_impute.tsv',
            'results/VIPER/liver_VIPER_impute.tsv',
            'results/ALRA/liver_ALRA_impute.tsv']

s_mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC',
          'SAVER', 'DCA', 'scIGANs', 'DrImpute', 'DeepImpute', 'scGNN', 'VIPER', 'ALRA']

labelPath = 'data/liver.stage.txt'
label = np.array(pd.read_csv(labelPath, header=None))

def generate_figures(dataPath, label, mtName, layouts, axs, type="paga", gene=None):
    raws, cols = layouts[0], layouts[1]
    for i in range(0, raws):
        for j in range(0, cols):
            print("generating {0}".format(mtName[cols * i + j]))
            data = pd.read_csv(dataPath[cols*i+j], sep = '\t', index_col = 0)
            data = pd.DataFrame(data.T, dtype = int)
            obs = pd.DataFrame(index=data.index)
            obs['label'] = label
            var_names = data.columns
            var = pd.DataFrame(index=var_names)
            adata = ad.AnnData(np.array(data), obs=obs, var=var)
            adata.X = adata.X.astype('float64')
            if type == "paga":
                sc.pp.filter_genes(adata, min_counts=1)
                sc.pp.recipe_zheng17(adata)
                sc.tl.pca(adata, svd_solver='arpack')
                sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                sc.tl.paga(adata, groups='label')
                sc.pl.paga(adata, threshold=0.07, cmap='gist_rainbow',show=False, fontsize=12, ax=axs[i, j], title=mtName[cols*i+j])
            elif type == "umap":
                sc.pp.filter_genes(adata, min_counts=1)
                sc.pp.recipe_zheng17(adata)
                sc.tl.pca(adata, svd_solver='arpack')
                sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                sc.tl.paga(adata, groups='label')
                sc.pl.paga(adata, threshold=0.07, cmap='gist_rainbow', show=False, fontsize=12,
                           title=mtName[cols * i + j])
                sc.tl.draw_graph(adata, init_pos='paga')
                sc.pl.draw_graph(adata, color='label', legend_loc='on data', ax=axs[i, j], show=False,
                                 title=mtName[cols * i + j])
            elif type == "violin":
                sc.pp.normalize_per_cell(adata) # the same as rescipe_zheng17
                sc.pp.log1p(adata)
                sc.pp.scale(adata)
                sc.pl.violin(adata, groupby='label', keys=[gene], ax=axs[i, j], show=False, stripplot=False, ylabel=gene)
                axs[i, j].title.set_text(mtName[cols*i+j])

# plot for paga
figa, axsa = plt.subplots(4, 3, figsize=(9,12),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (4, 3), axsa, type="paga")
figa.savefig("paper/liver/liver_paga1.pdf")

fige, axse = plt.subplots(4, 3, figsize=(9,12),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (4, 3), axse, type="umap")
fige.savefig("paper/liver/liver_umap1.pdf")