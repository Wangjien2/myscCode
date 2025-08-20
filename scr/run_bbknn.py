'''
Author: wangje
Date: 2025-08-11 10:12:52
'''
import scanpy as sc
import bbknn
import pandas as pd
import numpy as np
import logging
import os
import argparse
import sys
import numpy as np
from sklearn.decomposition import PCA
# pip install scikit-dimension
import skdim.id
# pip install kneed
from kneed import KneeLocator

## 读入数据
adata = sc.read_h5ad(args.input)

######## 数据预处理 ########
def preprocess_data(adata, target_sum, n_top_genes, max_value):
    """Preprocess AnnData: normalize, log-transform, select highly variable genes, regress out, and scale."""
    logging.info("Starting preprocessing")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    logging.info(f"Identified {sum(adata.var.highly_variable)} highly variable genes")
    
    # Save highest expressed genes plot
    sc.pl.highest_expr_genes(adata, n_top=20, save="_highest_expr_genes.png")
    
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ["percent.mt"])
    sc.pp.scale(adata, max_value=max_value)
    return adata

###### 获得最佳的Dim值 ######
def estimate_dimensionality(adata, var_threshold=0.95, use_pca=True, n_components=50):
    """
    用7种方法估计 AnnData 对象的内在维度，返回结果字典
    - adata: AnnData 对象，包含单细胞数据
    - var_threshold: PCA 方差解释阈值（如 0.95）
    - use_pca: 是否使用预计算的 PCA 结果（默认为 True）
    - n_components: 如果重新计算 PCA，使用的最大主成分数
    """
    results = {}

    # 选择输入数据
    if use_pca:
        # 检查是否已运行 PCA
        if 'X_pca' not in adata.obsm or 'pca' not in adata.uns:
            print("No PCA results found in adata. Running PCA...")
            sc.tl.pca(adata, n_comps=n_components)
        X = adata.obsm['X_pca']  # 使用 PCA 结果
        explained_variance = np.cumsum(adata.uns['pca']['variance_ratio'])
    else:
        X = adata.X  # 直接使用原始数据（可能需要预处理）
        # 如果使用原始数据，运行 PCA
        pca = PCA().fit(X)
        explained_variance = np.cumsum(pca.explained_variance_ratio_)
    var_dim = np.argmax(explained_variance >= var_threshold) + 1
    results['var_dim (PCA 95%)'] = var_dim
    mle = skdim.id.MLE().fit(X)
    results['mle_dim (MLE)'] = int(round(mle.dimension_))
    kn = KneeLocator(
        x=range(1, len(explained_variance) + 1),
        y=explained_variance,
        curve='convex', direction='decreasing',S=1.0  # 平滑参数，减少噪声干扰
    )
    elbow_dim = kn.knee if kn.knee is not None else var_dim
    results['elbow_dim (Elbow)'] = elbow_dim
    results['optimal_dim (Conservative)'] = min(
        var_dim, results['mle_dim (MLE)'], elbow_dim
    )
    twoNN = skdim.id.TwoNN().fit(X)
    results['twoNN_dim (TwoNN)'] = int(round(twoNN.dimension_))
    fishers = skdim.id.FisherS().fit(X)
    results['fisherS_dim (FisherS)'] = int(round(fishers.dimension_))
    ess = skdim.id.ESS().fit(X)
    results['ess_dim (ESS)'] = int(round(ess.dimension_))
    adata.uns['dimensionality_estimation'] = results
    return results

###### 获取最佳的n_neighbors值 ######
def estimate_n_neighbors(n_cells, complexity='medium', rare_cells=False, base='ln'):
    """
    估算单细胞分析中邻居图的 n_neighbors 值。

    参数:
        n_cells (int): 数据集中的细胞总数。
        complexity (str): 数据复杂性，'low'、'medium' 或 'high'。默认 'medium'。
        rare_cells (bool): 是否优先检测稀有细胞群体。默认 False。
        base (str): 对数基，'ln' 或 'log10'。默认 'ln'。

    返回:
        int: 推荐的 n_neighbors 值。
        dict: 包含推荐值、最小值、最大值的字典。
    """
    if not isinstance(n_cells, (int, float)) or n_cells < 1:
        raise ValueError("n_cells must be a positive number.")
    if base not in ['ln', 'log10']:
        raise ValueError("base must be 'ln' or 'log10'.")
    if complexity not in ['low', 'medium', 'high']:
        raise ValueError("complexity must be 'low', 'medium', or 'high'.")

    if base == 'ln':
        log_n = np.log(n_cells)
    else:
        log_n = np.log10(n_cells)

    base_factor = {'low': 2, 'medium': 3, 'high': 4}[complexity]
    n_neighbors = int(base_factor * log_n)

    if rare_cells:
        n_neighbors = int(n_neighbors * 0.7)

    min_neighbors = 10
    max_neighbors = 50
    if n_cells < 1000:
        min_neighbors = 5
        max_neighbors = 20
    elif n_cells > 50000:
        max_neighbors = 100

    n_neighbors = max(min_neighbors, min(n_neighbors, max_neighbors))

    result = {
        'recommended': n_neighbors,
        'min': min_neighbors,
        'max': max_neighbors
    }

    return n_neighbors, result

###### 批次校正和聚类 ######
logging.info("Running PCA")
sc.tl.pca(adata, svd_solver="arpack", mask_var="highly_variable")
sc.external.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs=22)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=22)
sc.tl.umap(adata, min_dist=min_dist, n_pcs=22)
sc.tl.leiden(adata)



adata = preprocess_data(adata, target_sum=1e4, n_top_genes=4000, max_value=10)