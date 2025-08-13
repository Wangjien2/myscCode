# harmony-pytorch地址：https://github.com/lilab-bcb/harmony-pytorch
# 安装：pip install harmony-pytorch
# 导入harmony-pytorch模块
import sys
import subprocess
import scanpy as sc
import scanpy.external as sce
from harmony import harmonize
import torch

# 检查GPU是否可用
if torch.cuda.is_available():
    # 如果可用，使用GPU
    device = torch.device("cuda")
    print("GPU is available.")
else:
    # 如果不可用，使用CPU
    device = torch.device("cpu")
    print("GPU is not available, using CPU instead.")

# 函数的用法
# def harmonize(
#     X: np.array,
#     df_obs: pd.DataFrame,
#     batch_key: Union[str, List[str]],
#     n_clusters: int = None,
#     max_iter_harmony: int = 10,
#     max_iter_clustering: int = 200,
#     tol_harmony: float = 1e-4,
#     tol_clustering: float = 1e-5,
#     ridge_lambda: float = 1.0,
#     sigma: float = 0.1,
#     block_proportion: float = 0.05,
#     init_centroids_method: str = "default",
#     theta: float = 2.0,
#     tau: int = 0,
#     random_state: int = 0,
#     use_gpu: bool = False,
#     n_jobs: int = -1,
#     verbose: bool = True,
# ) -> np.array:

# 读入数据
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


# harmony-pytorch
# harmony-pytorch模块的harmonize函数
# harmonize函数的参数说明：
# X: np.array，数据矩阵，每一行代表一个细胞，每一列代表一个基因
# df_obs: pd.DataFrame，观测数据，包含细胞的注释信息
# batch_key: Union[str, List[str]]，批处理键，用于指定批次信息
# n_clusters: int，聚类数目
# max_iter_harmony: int，Harmony算法的最大迭代次数
# max_iter_clustering: int，聚类算法的最大迭代次数
# tol_harmony: float，Harmony算法的收敛阈值
# tol_clustering: float，聚类算法的收敛阈值
# ridge_lambda: float，岭回归的参数
# sigma: float，高斯核函数的参数
# block_proportion: float，块的比例
# init_centroids_method: str，初始化质心的方法
# theta: float，Harmony算法的参数
# tau: int，Harmony算法的参数
# random_state: int，随机种子
# use_gpu: bool，是否使用GPU
# n_jobs: int，使用的线程数
# verbose: bool，是否打印详细信息

print("\n开始使用Harmony进行批次效应校正...")
# 使用harmony-pytorch进行批次效应校正
ho = harmonize(
    X=adata.obsm['X_pca'],
    df_obs=adata.obs, 
    batch_key='orig.ident', # 批次信息
    n_clusters=2,
    max_iter_harmony=10,
    max_iter_clustering=200,
    tol_harmony=1e-4,
    tol_clustering=1e-5,
    ridge_lambda=1.0,
    sigma=0.1,
    block_proportion=0.05,
    init_centroids_method="default",
    theta=2.0,
    tau=0,
    random_state=0,
    use_gpu=False,
    n_jobs=-1,
    verbose=True
)
# 将校正后的数据存储在obsm中
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
# 使用scanpy的neighbors和umap函数进行邻居图构建和UMAP降维
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
# sc.tl.tsne(adata, use_rep='X_pca_harmony')
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_harmony', use_weights=True)
# sc.tl.louvain(adata, resolution=0.5, key_added='louvain_harmony', use_weights=True)
# 将raw数据恢复到adata中
adata = adata.raw.to_adata()
print("\n恢复raw数据后的信息:")
print(f"raw数据包含 {adata.raw.n_obs} 个细胞和 {adata.raw.n_vars} 个基因")

# 保存结果
adata.write_h5ad("./scRNA_Harmony_harmony-pytorch.h5ad")
