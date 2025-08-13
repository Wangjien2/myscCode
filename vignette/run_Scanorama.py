'''
Author: wangje
Date: 2025-08-13 15:05:11
'''
# Scanorama官网地址:https://github.com/brianhie/scanorama
# Scanorama使用教程: https://github.com/brianhie/scanorama#full-tutorial
# 安装：pip install scanorama
import scanorama
import sys
import os
import scanpy as sc
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

# 执行Scanorama
sce.pp.scanorama_integrate(adata, key='orig.ident', basis='X_pca') 
print("\nScanorama后的数据信息:")
print(adata)

sc.pp.neighbors(adata, use_rep='X_scanorama')
sc.tl.umap(adata)
# sc.tl.tsne(adata, use_rep='X_scanorama')
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanorama',use_weights=True)
# sc.tl.louvain(adata, resolution=0.5, key_added='louvain_scanorama',use_weights=True)

# 将raw数据恢复到adata中
adata = adata.raw.to_adata()
print("\n恢复raw数据后的信息:")
print(f"raw数据包含 {adata.raw.n_obs} 个细胞和 {adata.raw.n_vars} 个基因")
# 保存数据
adata.write_h5ad("./scRNA_Scanorama.h5ad")
# 保存绘图数据
sc.pl.umap(adata, color=['leiden_scanorama'], save='Scanorama_leiden_scanorama.png')
sc.pl.umap(adata, color=['sample'], save='Scanorama_sample.png')
sc.pl.umap(adata, color=['orig.ident'], save='Scanorama_orig.ident.png')