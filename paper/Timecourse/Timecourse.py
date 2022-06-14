import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import anndata as ad
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
# dataPath = ['data/Timecourse.raw.tsv',
#             'data/Timecourse_SCDD_impute.tsv',
#             'data/Timecourse_Diffusion_impute.tsv',
#             'data/Timecourse_MAGIC_impute.tsv',
#             'data/Timecourse_SAVER_impute.tsv',
#             'data/Timecourse_DCA_impute.tsv',
#             'data/Timecourse_DeepImpute_impute.txt',
#             'data/Timecourse_SCRABBLE_impute.tsv',
#             'data/Timecourse_VIPER_impute.tsv',
#             'data/Timecourse_DrImpute_impute.tsv',
#             'data/Timecourse_scImpute_impute.tsv',
#             'data/Timecourse_scIGANs_impute.tsv']
# mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC',
#           'SAVER', 'DCA', 'DeepImpute', 'SCRABBLE',
#           'VIPER', 'DrImpute', 'scImpute', 'scIGANs']

s_dataPath = ['data/Timecourse.raw.tsv',
            'results/SCDD/Timecourse_SCDD_impute.tsv',
            'results/Diffusion/Timecourse_Diffusion_impute.tsv',
            'results/scImpute/Timecourse_scImpute_impute.tsv',
            'results/scIGANs/Timecourse_scIGANs_impute.tsv',
            'results/VIPER/Timecourse_VIPER_impute.tsv',
            'results/SAVER/Timecourse_SAVER_impute.tsv',
            'results/DCA/Timecourse_DCA_impute.tsv',
            'results/scVI/Timecourse_scVI_impute.tsv',
            'results/DrImpute/Timecourse_DrImpute_impute.tsv',
            'results/DeepImpute/Timecourse_DeepImpute_impute.tsv',
            'results/scGNN/Timecourse_scGNN_impute.tsv',
            'results/MAGIC/Timecourse_MAGIC_impute.tsv',
            'results/ALRA/Timecourse_ALRA_impute.tsv',
            'results/SCRABBLE/Timecourse_SCRABBLE_impute.tsv']

s_mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'scImpute', 'scIGANs', 'VIPER',
          'SAVER', 'DCA', 'scVI', 'DrImpute', 'DeepImpute', 'scGNN', 'MAGIC', 'ALRA', 'SCRABBLE']

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
                    sc.pl.paga(adata, threshold=None, cmap='gist_rainbow', show=False, fontsize=12,node_size_scale=0.7,
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


# # plot for paga
figa, axsa = plt.subplots(1, 3, figsize=(6,2),constrained_layout=True)
generate_figures(s_dataPath[:3], label, s_mtName[:3], (1, 3), axsa, type="paga")
figa.savefig("paper/Timecourse/Timecourse_paga31.svg")
#
figa, axsa = plt.subplots(3, 5, figsize=(10,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (3, 5), axsa, type="paga")
figa.savefig("paper/Timecourse/Timecourse_paga4.svg")
#
# fige, axse = plt.subplots(1, 3, figsize=(9,3),constrained_layout=True)
# generate_figures(s_dataPath[:3], label, s_mtName[:3], (1, 3), axse, type="umap")
# fige.savefig("paper/Timecourse/Timecourse_umap31.svg")
# #
# fige, axse = plt.subplots(3, 5, figsize=(15,9),constrained_layout=True)
# generate_figures(s_dataPath, label, s_mtName, (3, 5), axse, type="umap")
# fige.savefig("paper/Timecourse/Timecourse_umap4.svg")
# #
# # # plot for violins. marker genes: `NANOG`, `SOX2`, `CER1`
# fige, axse = plt.subplots(1, 3, figsize=(9,3),constrained_layout=True)
# generate_figures(s_dataPath[:3], label, s_mtName[:3], (1, 3), axse, type="violin", gene="NANOG")
# fige.savefig("paper/Timecourse/Timecourse_NANOG31.svg")
# # #
# fige, axse = plt.subplots(3, 5, figsize=(15,9),constrained_layout=True)
# generate_figures(s_dataPath, label, s_mtName, (3, 5), axse, type="violin", gene="NANOG")
# fige.savefig("paper/Timecourse/Timecourse_NANOG4.svg")
# # #
# fige, axse = plt.subplots(1, 3, figsize=(9,3),constrained_layout=True)
# generate_figures(s_dataPath[:3], label, s_mtName[:3], (1, 3), axse, type="violin", gene="SOX2")
# fige.savefig("paper/Timecourse/Timecourse_SOX231.svg")
# # #
# fige, axse = plt.subplots(3, 5, figsize=(15,9),constrained_layout=True)
# generate_figures(s_dataPath, label, s_mtName, (3, 5), axse, type="violin", gene="SOX2")
# fige.savefig("paper/Timecourse/Timecourse_SOX24.svg")
# #
# fige, axse = plt.subplots(1, 3, figsize=(9,3),constrained_layout=True)
# generate_figures(s_dataPath[:3], label, s_mtName[:3], (1, 3), axse, type="violin", gene="CER1")
# fige.savefig("paper/Timecourse/Timecourse_CER131.svg")
# #
# fige, axse = plt.subplots(3, 5, figsize=(15,9),constrained_layout=True)
# generate_figures(s_dataPath, label, s_mtName, (3, 5), axse, type="violin", gene="CER1")
# fige.savefig("paper/Timecourse/Timecourse_CER14.svg")

