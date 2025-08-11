if you just have a Seurat object in R. but you want run it in python by bbknn batch correct. you have to transfer seurat object int h5ad object by [sceasy](https://github.com/cellgeni/sceasy) or [schard](https://github.com/cellgeni/schard).

```python
set(adata.obs["orig.ident"].tolist()) adata = sc.read_h5ad("./scRNA_Harmony.h5ad")# import python moudle
import scanpy as sc
import pandas as pd
import bbknn
import os
import sys

# load data
adata = sc.read_h5ad("./scRNA_Harmony.h5ad")
adata
#AnnData object with n_obs × n_vars = 391885 × 33538
#    obs: 'orig.ident', 'percent.mt', 'percent.hb', 'percent.rb', 'scrublet_score', 'nCount_RNA', 'nFeature_RNA', #'Harmony_SNN_res.0.6'
#    var: 'name'
#    obsm: 'X_harmony', 'X_harmonypca'

# filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var_names_make_unique()
adata.obs_names_make_unique()

# save the raw data
adata.raw = adata.copy()

# Normalization and Feature Selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 188795 × 4457
#    obs: 'orig.ident', 'percent.mt', 'percent.hb', 'percent.rb', 'scrublet_score', 'nCount_RNA', 'nFeature_RNA', #'Harmony_SNN_res.0.6', 'n_genes'
#    var: 'name', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#    uns: 'log1p', 'hvg'
#    obsm: 'X_harmony', 'X_harmonypca'

# run PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# run BBKNN
## here I use orig.ident as batch column, you can choose you batch column
set(adata.obs["orig.ident"].tolist())
adata.obs[["orig.ident"]].value_counts()
bbknn.bbknn(adata, batch_key="orig.ident")

# cluster
sc.tl.louvain(adata)  # 结果保存在 adata.obs['louvain']
sc.tl.leiden(adata)   # 结果保存在 adata.obs['leiden']
# umap
sc.tl.umap(adata)
```
