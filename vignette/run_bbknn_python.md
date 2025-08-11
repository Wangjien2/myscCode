if you just have a Seurat object in R. but you want run it in python by bbknn batch correct. you have to transfer seurat object int h5ad object by [sceasy](https://github.com/cellgeni/sceasy) or [schard](https://github.com/cellgeni/schard).

```python
# 导入Python模块
import scanpy as sc
import pandas as pd
import bbknn
import os
import sys

# 加载数据
adata = sc.read_h5ad("./scRNA_Harmony.h5ad")
print("原始数据信息:")
print(adata)

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

# 运行BBKNN进行批次校正
print("\n批次信息分布:")
print(adata.obs["orig.ident"].value_counts())
bbknn.bbknn(adata, batch_key="orig.ident")

# 聚类分析
sc.tl.louvain(adata)
sc.tl.leiden(adata)

# 降维可视化（UMAP）
sc.tl.umap(adata)

# 确保输出目录存在
output_dir = "./processed_data"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 验证raw数据是否存在
if adata.raw is not None:
    print("\n准备保存数据：raw数据已确认存在")
else:
    print("\n警告：保存前未检测到raw数据！")

# 保存完整的AnnData对象（包含raw数据）
output_path = os.path.join(output_dir, "scRNA_processed.h5ad")
adata.write(output_path)
print(f"处理完成的数据集已保存至: {output_path}")

# 加载数据验证raw是否被正确保存
loaded_adata = sc.read_h5ad(output_path)
if loaded_adata.raw is not None:
    print("\n保存验证：raw数据成功保留")
    print(f"加载后raw数据包含 {loaded_adata.raw.n_obs} 个细胞和 {loaded_adata.raw.n_vars} 个基因")
else:
    print("\n保存验证：警告！加载后未发现raw数据")

# 保存细胞聚类信息
cluster_info = loaded_adata.obs[['orig.ident', 'louvain', 'leiden', 'percent.mt']]
cluster_info.to_csv(os.path.join(output_dir, "cell_cluster_info.csv"))
print(f"细胞聚类信息已保存至: {os.path.join(output_dir, 'cell_cluster_info.csv')}")

```

After batch correction and clustering, we can annotate the cell type using some automatic annotation software.
