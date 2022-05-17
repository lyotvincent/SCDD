import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import anndata as ad
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
dataPath = ['data/Timecourse.raw.tsv', 
            'data/Timecourse_SCDD_impute.tsv',
            'data/Timecourse_Diffusion_impute.tsv',
            'data/Timecourse_MAGIC_impute.tsv',
            'data/Timecourse_SAVER_impute.tsv',
            'data/Timecourse_DCA_impute.tsv',
            'data/Timecourse_DeepImpute_impute.txt',
            'data/Timecourse_SCRABBLE_impute.tsv',
            'data/Timecourse_VIPER_impute.tsv',
            'data/Timecourse_DrImpute_impute.tsv',
            'data/Timecourse_scImpute_impute.tsv', 
            'data/Timecourse_scIGANs_impute.tsv']
mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC',
          'SAVER', 'DCA', 'DeepImpute', 'SCRABBLE',
          'VIPER', 'DrImpute', 'scImpute', 'scIGANs']

s_dataPath = ['data/Timecourse.raw.tsv', 
            'results/SCDD/Timecourse_SCDD3_impute.tsv',
            'results/Diffusion/Timecourse_Diffusion_impute.tsv',
            'results/MAGIC/Timecourse_MAGIC_impute.tsv',
            'results/SAVER/Timecourse_SAVER_impute.tsv',
            'results/DCA/Timecourse_DCA_impute.tsv']

s_mtName = ['Raw', 'SCDD', 'SCDD(Diffusion)', 'MAGIC',
          'SAVER', 'DCA']

labelPath = 'data/Timecourse.label.txt'
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
            sc.pp.filter_genes(adata, min_counts=1)
            sc.pp.recipe_zheng17(adata)
            sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
            sc.tl.paga(adata, groups='label')
            if type == "paga":
                sc.pl.paga(adata, threshold=None, cmap='gist_rainbow',show=False, fontsize=12, ax=axs[i, j], title=mtName[cols*i+j])
            elif type == "violin":
                sc.pl.violin(adata, groupby='label', keys=[gene], ax=axs[i, j], show=False, stripplot=False, ylabel=gene)
                axs[i, j].title.set_text(mtName[cols*i+j])

# plot for paga
figa, axsa = plt.subplots(2, 3, figsize=(9,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (2, 3), axsa, type="paga")

# plot for violins. marker genes: `NANOG`, `SOX2`, `CER1`
figb, axsb = plt.subplots(2, 3, figsize=(9,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (2, 3), axsb, type="violin", gene="NANOG")

figc, axsc = plt.subplots(2, 3, figsize=(9,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (2, 3), axsc, type="violin", gene="SOX2")

figd, axsd = plt.subplots(2, 3, figsize=(9,6),constrained_layout=True)
generate_figures(s_dataPath, label, s_mtName, (2, 3), axsd, type="violin", gene="CER1")

# save figures
figa.savefig("paper/Timecourse/Timecourse_paga.pdf")
figb.savefig("paper/Timecourse/Timecourse_nanog.pdf")
figc.savefig("paper/Timecourse/Timecourse_sox2.pdf")
figd.savefig("paper/Timecourse/Timecourse_cer1.pdf")
