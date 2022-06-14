import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import anndata as ad
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


mask_dataPath = ['data/Timecourse.raw.tsv',
              'data/mask0.1_Timecourse.raw.tsv',
              'data/mask0.2_Timecourse.raw.tsv',
              'data/mask0.3_Timecourse.raw.tsv',
              'data/mask0.4_Timecourse.raw.tsv',
              'results/SCDD/Timecourse_SCDD_impute.tsv',
              'results/SCDD/mask0.1Timecourse_SCDD1_impute.tsv',
              'results/SCDD/mask0.2Timecourse_SCDD1_impute.tsv',
              'results/SCDD/mask0.3Timecourse_SCDD1_impute.tsv',
              'results/SCDD/mask0.4Timecourse_SCDD1_impute.tsv']

mask_mtName = ['Raw', 'mask 10%', 'mask 20%', 'mask 30%', 'mask 40%',
               'SCDD', 'SCDD 10%', 'SCDD 20%', 'SCDD 30%', 'SCDD 40%']

labelPath = 'data/Timecourse.label.txt'
label = np.array(pd.read_csv(labelPath, header=None))

def generate_figures(dataPath, label, mtName, layouts, axs, type="paga", gene=None):
    raws, cols = layouts[0], layouts[1]
    if raws == 1:
        for i in range(0, cols):
            if mtName[i] == 'scImpute':
                data = pd.read_csv(dataPath[i], sep=' ')
            else:
                data = pd.read_csv(dataPath[i], sep='\t', index_col=0)
            print("generating {0}".format(mtName[i]))
            data = pd.DataFrame(data.T, dtype=int)
            obs = pd.DataFrame(index=data.index)
            obs[''] = label
            var_names = data.columns
            var = pd.DataFrame(index=var_names)
            adata = ad.AnnData(np.array(data), obs=obs, var=var)
            adata.X = adata.X.astype('float64')
            if type == "paga":
                sc.pp.filter_genes(adata, min_counts=1)
                sc.pp.recipe_zheng17(adata)
                sc.tl.pca(adata, svd_solver='arpack')
                sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                sc.tl.paga(adata, groups='')
                sc.pl.paga(adata, threshold=None, cmap='gist_rainbow', show=False, fontsize=12, ax=axs[i],
                           title=mtName[i])
            elif type == "umap":
                sc.pp.filter_genes(adata, min_counts=1)
                sc.pp.recipe_zheng17(adata)
                sc.tl.pca(adata, svd_solver='arpack')
                sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                sc.tl.paga(adata, groups='')
                sc.pl.paga(adata, threshold=None, cmap='gist_rainbow', show=False, fontsize=12,
                           title=mtName[i])
                sc.tl.draw_graph(adata, init_pos='paga')
                sc.pl.draw_graph(adata, color='', legend_loc='on data', ax=axs[i], show=False,
                                 title=mtName[i])
            elif type == "violin":
                sc.pp.normalize_per_cell(adata)  # the same as rescipe_zheng17
                sc.pp.log1p(adata)
                sc.pp.scale(adata)
                sc.pl.violin(adata, groupby='', keys=[gene], ax=axs[i], show=False, stripplot=False,
                             ylabel=gene)
                axs[i].title.set_text(mtName[i])
    else:
        for i in range(0, raws):
            for j in range(0, cols):
                print("generating {0}".format(mtName[cols * i + j]))
                if mtName[cols * i + j] == 'scImpute':
                    data = pd.read_csv(dataPath[cols*i+j], sep=' ')
                else:
                    data = pd.read_csv(dataPath[cols*i+j], sep = '\t', index_col = 0)
                data = pd.DataFrame(data.T, dtype = int)
                obs = pd.DataFrame(index=data.index)
                obs[''] = label
                var_names = data.columns
                var = pd.DataFrame(index=var_names)
                adata = ad.AnnData(np.array(data), obs=obs, var=var)
                adata.X = adata.X.astype('float64')
                if type == "paga":
                    sc.pp.filter_genes(adata, min_counts=1)
                    sc.pp.recipe_zheng17(adata)
                    sc.tl.pca(adata, svd_solver='arpack')
                    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                    sc.tl.paga(adata, groups='')
                    sc.pl.paga(adata, threshold=None, cmap='gist_rainbow',show=False, fontsize=12, ax=axs[i, j], title=mtName[cols*i+j])
                elif type == "umap":
                    sc.pp.filter_genes(adata, min_counts=1)
                    sc.pp.recipe_zheng17(adata)
                    sc.tl.pca(adata, svd_solver='arpack')
                    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
                    sc.tl.paga(adata, groups='')
                    sc.pl.paga(adata, threshold=None, cmap='gist_rainbow', show=False, fontsize=12,
                               title=mtName[cols * i + j])
                    sc.tl.draw_graph(adata, init_pos='paga')
                    sc.pl.draw_graph(adata, color='', legend_loc='on data', ax=axs[i, j], show=False,
                                     title=mtName[cols * i + j])
                elif type == "violin":
                    sc.pp.normalize_per_cell(adata) # the same as rescipe_zheng17
                    sc.pp.log1p(adata)
                    sc.pp.scale(adata)
                    sc.pl.violin(adata, groupby='', keys=[gene], ax=axs[i, j], show=False, stripplot=False, ylabel=gene)
                    axs[i, j].title.set_text(mtName[cols*i+j])

figa, axsa = plt.subplots(2, 5, figsize=(15,6),constrained_layout=True)
generate_figures(mask_dataPath, label, mask_mtName, (2, 5), axsa, type="paga")
figa.savefig("paper/Timecourse/Timecourse_mask_paga.svg")

figa, axsa = plt.subplots(2, 5, figsize=(15,6),constrained_layout=True)
generate_figures(mask_dataPath, label, mask_mtName, (2, 5), axsa, type="umap")
figa.savefig("paper/Timecourse/Timecourse_mask_umap.svg")

figa, axsa = plt.subplots(2, 5, figsize=(15,6),constrained_layout=True)
generate_figures(mask_dataPath, label, mask_mtName, (2, 5), axsa, type="violin", gene='NANOG')
figa.savefig("paper/Timecourse/Timecourse_mask_NANOG.svg")

figa, axsa = plt.subplots(2, 5, figsize=(15,6),constrained_layout=True)
generate_figures(mask_dataPath, label, mask_mtName, (2, 5), axsa, type="violin", gene='SOX2')
figa.savefig("paper/Timecourse/Timecourse_mask_SOX2.svg")

figa, axsa = plt.subplots(2, 5, figsize=(15,6),constrained_layout=True)
generate_figures(mask_dataPath, label, mask_mtName, (2, 5), axsa, type="violin", gene='CER1')
figa.savefig("paper/Timecourse/Timecourse_mask_CER1.svg")