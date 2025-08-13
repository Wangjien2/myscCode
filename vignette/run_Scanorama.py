'''
Author: wangje
Date: 2025-08-13 15:05:11
'''# Scanorama官网地址:https://github.com/brianhie/scanorama
# Scanorama使用教程: https://github.com/brianhie/scanorama#full-tutorial
# 安装：pip install scanorama
import scanorama
import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd
import scanpy.external as sce

# 读取数据
adata = sc.read_h5ad("./scRNA_Harmony.h5ad")
# 数据过滤
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var_names_make_unique()
adata.obs_names_make_unique()

# 保存原始数据到raw属性
adata.raw = adata.copy()
print("\n设置raw数据后的信息:")
print(f"raw数据包含 {adata.raw.n_obs} 个细胞和 {adata.raw.n_vars} 个基因")

# 标准化和特征选择
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
print("\n筛选高变基因后的数据信息:")
print(adata)

# 执行PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
print("\nPCA后的数据信息:")
print(adata)

# 执行Scanorama整合
# 保存原始基因名称
gene_names = adata.var_names.tolist()

# 按批次拆分数据集并保留相同基因
batches = adata.obs['orig.ident'].unique()
datasets = []
for batch in batches:
    # 获取当前批次数据并确保基因顺序一致
    batch_data = adata[adata.obs['orig.ident'] == batch, :]
    # 确保所有批次使用相同的基因集合
    datasets.append(batch_data[:, gene_names].X)

# 执行Scanorama整合
integrated, genes = scanorama.integrate(
    datasets,
    gene_names,
    dimred=50,
    verbose=True
)

# 将整合结果合并回AnnData对象
adata.obsm['X_scanorama'] = np.concatenate(integrated)
print("\nScanorama后的数据信息:")
print(adata)

sc.pp.neighbors(adata, use_rep='X_scanorama')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanorama', use_weights=True)

# 保存数据（不恢复raw数据以保留整合结果）
adata.write_h5ad("./scRNA_Scanorama.h5ad")