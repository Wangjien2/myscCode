# Function: RunBBKNNClustering
# Description: Perform BBKNN batch correction and clustering on a Seurat object
# Parameters:
#   srt: Seurat object containing single-cell data
#   batch: Column name in metadata specifying batch information
#   assay: Assay to use (default: DefaultAssay(srt))
#   linear_reduction: Linear reduction method (default: "pca")
#   linear_reduction_dims: Number of dimensions for linear reduction (default: 50)
#   linear_reduction_dims_use: Specific dimensions to use for BBKNN (default: NULL, auto-estimated)
#   cluster_algorithm: Clustering algorithm ("louvain", "slm", or "leiden") (default: "louvain")
#   cluster_resolution: Resolution for clustering (default: 0.6)
#   bbknn_params: Additional parameters for BBKNN (default: list())
#   seed: Random seed for reproducibility (default: 11)
# Returns: Seurat object with BBKNN graph and clustering results

RunBBKNNClustering <- function(srt, batch, assay = NULL, 
                              linear_reduction = "pca", linear_reduction_dims = 50, 
                              linear_reduction_dims_use = NULL, 
                              cluster_algorithm = "louvain", cluster_resolution = 0.6, 
                              bbknn_params = list(), seed = 11) {
  library(Seurat)
  library(reticulate)
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Validate inputs
  if (!inherits(srt, "Seurat")) {
    stop("Input 'srt' must be a Seurat object.")
  }
  if (!batch %in% colnames(srt@meta.data)) {
    stop("Batch column '", batch, "' not found in Seurat object metadata.")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (!assay %in% Assays(srt)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt))
  if (!linear_reduction %in% reduc_test) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    if (!requireNamespace("leidenAlg", quietly = TRUE)) {
      stop("Package 'leidenAlg' is required for Leiden clustering. Please install it.")
    }
  }
  
  # Check Python BBKNN availability
  if (!py_module_available("bbknn")) {
    stop("Python module 'bbknn' is not available. Please install it using pip.")
  }
  bbknn <- import("bbknn")
  
  # Ensure variable features are available
  if (length(VariableFeatures(srt, assay = assay)) == 0) {
    message("No variable features found. Running FindVariableFeatures...")
    srt <- FindVariableFeatures(srt, assay = assay, nfeatures = 2000, verbose = FALSE)
  }
  
  # Scale data if not already scaled
  if (!all(VariableFeatures(srt, assay = assay) %in% rownames(GetAssayData(srt, slot = "scale.data", assay = assay)))) {
    message("Scaling data...")
    srt <- ScaleData(srt, assay = assay, features = VariableFeatures(srt, assay = assay), verbose = FALSE)
  }
  
  # Perform linear dimension reduction
  message("Performing linear dimension reduction (", linear_reduction, ")...")
  srt <- RunPCA(srt, assay = assay, features = VariableFeatures(srt, assay = assay), 
                npcs = linear_reduction_dims, reduction.name = "BBKNNpca", 
                verbose = FALSE, seed.use = seed)
  
  # Estimate dimensions if not provided
  if (is.null(linear_reduction_dims_use)) {
    dim_est <- tryCatch({
      min(
        intrinsicDimension::maxLikGlobalDimEst(
          data = Embeddings(srt, reduction = "BBKNNpca"), k = 20
        )[["dim.est"]],
        linear_reduction_dims
      )
    }, error = function(e) {
      message("Cannot estimate intrinsic dimensions. Using default range...")
      return(min(linear_reduction_dims, 30))
    })
    linear_reduction_dims_use <- seq_len(max(10, ceiling(dim_est)))
  }
  
  # Run BBKNN
  message("Running BBKNN integration...")
  emb <- Embeddings(srt, reduction = "BBKNNpca")[, linear_reduction_dims_use, drop = FALSE]
  params <- list(
    pca = emb,
    batch_list = srt[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke(.fn = bbknn$matrix$bbknn, .args = params)
  n.neighbors <- bem[[3]]$n_neighbors
  
  # Store BBKNN graph
  bbknn_graph <- as.sparse(bem[[2]][1:nrow(bem[[2]]), , drop = FALSE])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- assay
  srt@graphs[["BBKNN"]] <- bbknn_graph
  
  # Store BBKNN distance matrix
  bbknn_dist <- t(as.sparse(bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]))
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- assay
  srt@graphs[["BBKNN_dist"]] <- bbknn_dist
  
  # Store BBKNN neighbors
  val <- split(bbknn_dist@x, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  pos <- split(bbknn_dist@i + 1, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  idx <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = FALSE), n.neighbors - 1)]
    length(out) <- n.neighbors - 1
    return(out)
  }, x = val, y = pos))
  idx[is.na(idx)] <- sample(seq_len(nrow(idx)), size = sum(is.na(idx)), replace = TRUE)
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = FALSE), n.neighbors - 1)]
    length(out) <- n.neighbors - 1
    out[is.na(out)] <- 0
    return(out)
  }, x = val, y = val))
  dist <- cbind(0, dist)
  srt[["BBKNN_neighbors"]] <- new(Class = "Neighbor", nn.idx = idx, nn.dist = dist, 
                                  alg.info = list(), cell.names = rownames(emb))
  
  # Perform clustering
  message("Performing clustering (", cluster_algorithm, ")...")
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
                                   louvain = 1,
                                   slm = 3,
                                   leiden = 4)
  srt <- tryCatch({
    srt <- FindClusters(srt, graph.name = "BBKNN", resolution = cluster_resolution, 
                        algorithm = cluster_algorithm_index, method = "igraph", verbose = FALSE)
    srt[["BBKNNclusters"]] <- Idents(srt)
    srt
  }, error = function(e) {
    message("Error during clustering: ", e$message, ". Skipping clustering...")
    return(srt)
  })
  
  # Set default assay
  DefaultAssay(srt) <- assay
  
  return(srt)
}