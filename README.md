# SCDD 

## Abstract

Single-cell sequencing technologies are widely used to discover the evolutionary relationships and the differences in cells. Since dropout events may frustrate the analysis, many imputation approaches for single-cell RNA-seq data have appeared in previous attempts. However, previous imputation attempts usually suffer from the over-smooth problem, which may bring limited improvement or negative effect for the downstream analysis of single-cell RNA-seq data. In order to solve this difficulty, we propose a novel two-stage diffusion-denoising method called SCDD for single-cell RNA-seq imputation in this paper. We introduce the diffusion, i.e., a direct imputation strategy using the expression of similar cells for potential dropout sites, to perform the initial imputation at first. After the diffusion, a joint model integrated with graph convolutional neural network and contractive autoencoder is developed to generate superposition states of similar cells, from which we restore the original states and remove the noise introduced by the diffusion. The final experimental results indicate that, SCDD could effectively suppress the over-smooth problem and remarkably improve the effect of single-cell RNA-seq downstream analysis, including clustering and trajectory analysis.

## Tutorial

### Installation

**For Python 3.8:**

The Python packages used in `SCDD` are in `requirements_SCDD.txt`.

```shell
pip install -r requirements_SCDD.txt
```

If your are going to generate the figures in this paper, you need to installed the Python packages in `requirements.txt`.

```shell
pip install -r requirements.txt
```

**For R 4.1.0:**

Run this code in R Console for `SCDD` requirements if you need to select the param `neighbor_method` as `SC3`:

```R
install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment","SC3"))
```

If your are going to generate the figures in this paper, run this code first:

```R
install.packages(c("ggplot2", "ggpubr", "ROCR", "aricode", "BiocManager"))
BiocManager::install(c("SingleCellExperiment", "SC3", "Seurat", "SeuratObject", "SAVER"))
```

### Usage

Run ***run_SCDD.py*** to begin imputation:

```
python run_SCDD.py
```

The meaning of optional argument in ***run_SCDD.py*** is listed below:
  **-h, --help**: show this help message and exit
  **-n, --name**:  the name of experiments or dataset to identify the imputed results
  **-r, --raw**: the raw data path, except format: `.tsv` and `.h5ad`, make sure that your delimiter in .tsv file is  and your raw data is .X in .h5ad file.
  **-tr, --trans**: if raw data layout is: columns are cells and rows are genes, set this param `True`, else `False`.
  **-i, --id:** the id number to identify different results from the same methods, default 0.
  **-f, --format**: the format of input raw data file, support `tsv` and `h5ad`, default None, this will automatically identify the format by the suffix of the raw file.
  **-me, --max-epoch**: the max epoch of Denoising, default 500.
  **-nb, --neighbors**: the neighbors of each cell, default 20.
  **-t, --threshold**: the threshold of whether the point is a drop-out point, the drop-out rate higher than this value will be treated as a drop-out point, default 0.2.
  **-nm, --neighbor-method**: the method to generate concensus matrix and neighbors relations, default `SC3`. if the number of cells is more than 10000, it will automatically turn to `SNN`.
  **-b, --batch-size**: the batch_size of each input in the nerual network when denoising, default 5000, which means if cell number are less than 5000, the total batch will be 1.

### Examples:

Using SCDD to impute `Li` data, which format is `tsv`:

```shell
python run_SCDD.py -n Li -r data/Li.raw.tsv -i 0 -f tsv -nm SC3 -me 2000 -nb 20 -t 0.2 
```

And the results will be in: `results/SCDD/Li_SCDD_impute.tsv` .

Using SCDD to impute `Bladder` data, which format is `h5ad`:

```shell
python run_SCDD.py -n Bladder -r data/TS_Bladder.h5ad -i 0 -nm SNN -me 400 -nb 20 -t 0.2
```

And the results will be in: `results/SCDD/Bladder_SCDD_impute.h5ad`.

### Figures

Figures are in the fold `paper`. If you are going to generate figures yourself, You can use `genereate` module in `main.py`:

```python
from generates import *
Generate_Cellcycle()
```
