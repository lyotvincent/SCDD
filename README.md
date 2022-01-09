# SCDD 
Single-cell sequencing technologies are widely used to discover the evolutionary relationships and the differences in cells. Since dropout events may frustrate the analysis, many imputation approaches for single-cell RNA-seq data have appeared in previous attempts. However, previous imputation attempts usually suffer from the over-smooth problem, which may bring limited improvement or negative effect for the downstream analysis of single-cell RNA-seq data. In order to solve this difficulty, we propose a novel two-stage diffusion-denoising method called SCDD for single-cell RNA-seq imputation in this paper. We introduce the diffusion, i.e., a direct imputation strategy using the expression of similar cells for potential dropout sites, to perform the initial imputation at first. After the diffusion, a joint model integrated with graph convolutional neural network and contractive autoencoder is developed to generate superposition states of similar cells, from which we restore the original states and remove the noise introduced by the diffusion. The final experimental results indicate that, SCDD could effectively suppress the over-smooth problem and remarkably improve the effect of single-cell RNA-seq downstream analysis, including clustering and trajectory analysis.

## Tutorial
Now, we have four default experiments:Cellcycle, Timecourse, Li and Bone. You can use SCDD like this:
```python
from SCDD_api import *
e = SCDD(name="Timecourse")
e.run()
```
and the Diffusion and Denoising results will be saved in the fold `results`. Or you can use the raw data as the input:
```python
from SCDD_api import *
e = SCDD(raw="data/Timecourse.raw.tsv")
e.run()
```


