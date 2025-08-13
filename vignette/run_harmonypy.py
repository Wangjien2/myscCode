import os
import sys
import subprocess
import scanpy as sc
import scanpy.external as sce

# 检查并安装harmonypy模块
try:
    import harmonypy as hm
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy"])
    import harmonypy as hm

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

# 运行Harmony进行批次效应校正
# 使用harmonypy进行批次效应校正,然后将校正后的数据存储在obsm中
# print("\n开始使用Harmony进行批次效应校正...")
# ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'orig.ident', max_iter_harmony=30)
# adata.obsm['X_pca_harmony'] = ho.Z_corr.T  # 将校正后的数据存储在obsm中

# 使用scanpy的harmony_integrate函数进行批次效应校正
print("\n开始使用Harmony进行批次效应校正...")
sce.pp.harmony_integrate(adata, key='orig.ident', max_iter_harmony=30)
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
# sc.tl.tsne(adata, use_rep='X_pca_harmony')
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_harmony',use_weights=True)
# sc.tl.louvain(adata, resolution=0.5, key_added='louvain_harmony',use_weights=True)

# 将raw数据恢复到adata中
adata = adata.raw.to_adata()
print("\n恢复raw数据后的信息:")
print(f"raw数据包含 {adata.raw.n_obs} 个细胞和 {adata.raw.n_vars} 个基因")

# 保存结果
adata.write_h5ad("./scRNA_Harmony_corrected.h5ad")
