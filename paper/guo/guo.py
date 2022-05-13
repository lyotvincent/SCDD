import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import anndata as ad
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
s_dataPath = ['data/guo.raw.txt', 
            'results/SCDD/guo_SCDD1_impute.tsv',
            'results/Diffusion/guo_Diffusion_impute.tsv',
            'results/MAGIC/guo_MAGIC_impute.tsv',
            'results/SAVER/guo_SAVER_impute.tsv',
            'results/DCA/guo_DCA_impute.tsv']

s_mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC',
          'SAVER', 'DCA']

labelPath = 'data/guo.label.txt'
label = np.array(pd.read_csv(labelPath, header=None))

def generate_figures(dataPath, label, mtName, layouts, axs, type="paga", gene=None):
    raws, cols = layouts[0], layouts[1]
    for i in range(0, raws):
        for j in range(0, cols):
            data = pd.read_csv(dataPath[cols*i+j], sep = '\t', index_col = 0)
            data = pd.DataFrame(data.T, dtype = int)
            obs = pd.DataFrame(index=data.index)
            obs['label'] = label
            var_names = data.columns
            var = pd.DataFrame(index=var_names)
            adata = ad.AnnData(np.array(data), obs=obs, var=var)
            adata.X = adata.X.astype('float64')
            sc.pp.recipe_zheng17(adata, n_top_genes=43)
            sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
            sc.tl.paga(adata, groups='label')
            if type == "paga":
                sc.pl.paga(adata, threshold=None, cmap='gist_rainbow',show=False, fontsize=12, title=mtName[cols*i+j])
                sc.tl.draw_graph(adata, init_pos='paga')
                sc.pl.draw_graph(adata, color='label', legend_loc='on data', ax=axs[i, j], show=False,
                                 title=mtName[cols*i+j])
            elif type == "violin":
                sc.pl.violin(adata, groupby='label', keys=[gene], ax=axs[i, j], show=False, stripplot=False, ylabel=gene)
                axs[i, j].title.set_text(mtName[cols*i+j])


# plot for paga
figa, axsa = plt.subplots(2, 3, figsize=(9,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (2, 3), axsa, type="paga")
figa.savefig("paper/guo/guo_umap.pdf")