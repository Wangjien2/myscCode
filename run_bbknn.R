#' Advanced Batch-aware BBKNN Clustering for Single-Cell RNA-seq Data
#' Fixed version: Resolves BBKNN parameter name mismatch error
#'
#' @param srt Seurat object containing single-cell transcriptomic data.
#' @param batch Column name in metadata with batch information.
#' @param assay Assay to use (default: active assay).
#' @param normalization.method Normalization method (default: "LogNormalize").
#' @param scale.factor Scale factor for normalization (default: 10000).
#' @param var.features.n Number of variable features per batch (default: 2000).
#' @param var.features.method Variable feature selection method (default: "vst").
#' @param recompute.normalization Force recompute normalization (default: FALSE).
#' @param recompute.var.features Force recompute variable features (default: FALSE).
#' @param linear.reduction Linear reduction method ("pca", "svd", "ica") (default: "pca").
#' @param linear.reduction.dims Total linear dimensions (default: 50).
#' @param linear.reduction.dims.use Dimensions to use (auto-estimated if NULL).
#' @param dim.est.k Neighbors for dimension estimation (default: 20).
#' @param bbknn.n_neighbors Total neighbors per cell (default: 100)
#' @param bbknn.n_neighbors_within_batch Same-batch neighbors (default: 3)
#' @param bbknn.params Additional BBKNN parameters (default: list())
#' @param cluster.algorithm Clustering algorithm ("louvain", "slm", "leiden") (default: "louvain").
#' @param cluster.resolution Clustering resolution (default: 0.6).
#' @param nonlin.reduction Non-linear reduction ("umap", "tsne", NULL) (default: "umap").
#' @param nonlin.dims Non-linear dimensions (default: 2).
#' @param umap.min.dist UMAP minimum distance (default: 0.3).
#' @param umap.n.neighbors UMAP neighbors (default: 30).
#' @param seed Random seed (default: 11).
#' @param verbose Show messages (default: TRUE).
#'
#' @return Seurat object with integrated results
#' @export
#' @importFrom Seurat DefaultAssay Assays VariableFeatures FindVariableFeatures ScaleData
#' @importFrom Seurat RunPCA RunSVD RunICA Reductions FindClusters Idents RunUMAP RunTSNE
#' @importFrom Seurat NormalizeData subset MergeSeurat Embeddings
#' @importFrom reticulate import py_module_available py_get_attr
#' @importFrom purrr exec
#' @importFrom methods as new
#' @importFrom utils modifyList
#' @importFrom Matrix sparseMatrix

RunBBKNNClustering <- function(
  srt, batch, assay = NULL,
  # Preprocessing parameters
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  var.features.n = 2000,
  var.features.method = "vst",
  recompute.normalization = FALSE,
  recompute.var.features = FALSE,
  # Dimensionality reduction parameters
  linear.reduction = "pca",
  linear.reduction.dims = 50,
  linear.reduction.dims.use = NULL,
  dim.est.k = 20,
  # BBKNN parameters
  bbknn.n_neighbors = 100,
  bbknn.n_neighbors_within_batch = 3,
  bbknn.params = list(),
  # Clustering parameters
  cluster.algorithm = "louvain",
  cluster.resolution = 0.6,
  # Visualization parameters
  nonlin.reduction = "umap",
  nonlin.dims = 2,
  umap.min.dist = 0.3,
  umap.n.neighbors = 30,
  # General parameters
  seed = 11,
  verbose = TRUE
) {
  # Set seed
  set.seed(seed)
  
  # --------------------------
  # Validate dependencies
  # --------------------------
  if (!py_module_available("bbknn")) {
    stop("Python module 'bbknn' is required. Install with: pip install bbknn")
  }
  
  bbknn <- import("bbknn")
  bbknn_version <- tryCatch({
    as.numeric_version(py_get_attr(bbknn, "__version__"))
  }, error = function(e) {
    if (verbose) message("Could not retrieve BBKNN version.")
    numeric_version("0.0.0")
  })
  
  if (bbknn_version < "1.5.0") {
    warning("BBKNN version ", bbknn_version, " detected. Recommend upgrading to 1.5+ for parameter compatibility: pip install --upgrade bbknn")
  }
  
  required.pkgs <- c("Seurat", "reticulate", "purrr", "methods", "Matrix")
  missing.pkgs <- required.pkgs[!sapply(required.pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing.pkgs) > 0) {
    stop("Missing required packages: ", paste(missing.pkgs, collapse = ", "))
  }
  
  if (cluster.algorithm == "leiden" && !requireNamespace("leiden", quietly = TRUE)) {
    stop("Package 'leiden' is required for Leiden clustering. Install with: install.packages('leiden')")
  }
  
  if (is.null(linear.reduction.dims.use) && !requireNamespace("intrinsicDimension", quietly = TRUE)) {
    stop("Package 'intrinsicDimension' required for auto-dimension estimation. Install with: install.packages('intrinsicDimension')")
  }
  
  # --------------------------
  # Validate parameters
  # --------------------------
  # Normalization and variable features
  valid.norm.methods <- c("LogNormalize", "SCTransform")
  if (!normalization.method %in% valid.norm.methods) {
    stop("Invalid 'normalization.method'. Must be one of: ", paste(valid.norm.methods, collapse = ", "))
  }
  
  valid.var.methods <- c("vst", "mean.var.plot", "dispersion")
  if (!var.features.method %in% valid.var.methods) {
    stop("Invalid 'var.features.method'. Must be one of: ", paste(valid.var.methods, collapse = ", "))
  }
  
  # BBKNN parameters
  if (!is.numeric(bbknn.n_neighbors) || bbknn.n_neighbors < 1 || bbknn.n_neighbors %% 1 != 0) {
    stop("'bbknn.n_neighbors' must be a positive integer")
  }
  if (!is.numeric(bbknn.n_neighbors_within_batch) || bbknn.n_neighbors_within_batch < 1 || bbknn.n_neighbors_within_batch %% 1 != 0) {
    stop("'bbknn.n_neighbors_within_batch' must be a positive integer")
  }
  if (bbknn.n_neighbors_within_batch > bbknn.n_neighbors) {
    stop("'bbknn.n_neighbors_within_batch' cannot exceed 'bbknn.n_neighbors'")
  }
  
  # Validate bbknn.params
  deprecated.params <- c("neighbors", "neighbors_within_batch")  # Old camelCase or wrong names
  if (any(deprecated.params %in% names(bbknn.params))) {
    stop("Deprecated BBKNN parameters detected in 'bbknn.params'. Use snake_case: 'n_neighbors', 'n_neighbors_within_batch'")
  }
  
  valid.bbknn.params <- c("metric", "n_pcs", "trim", "annoy_n_trees", "n_trees", "use_faiss", "set_op_mix_ratio")
  invalid.params <- setdiff(names(bbknn.params), valid.bbknn.params)
  if (length(invalid.params) > 0) {
    warning("Invalid or unsupported BBKNN parameters provided: ", paste(invalid.params, collapse = ", "), ". They may cause errors.")
  }
  
  # Clustering
  if (!is.numeric(cluster.resolution) || cluster.resolution <= 0) {
    stop("'cluster.resolution' must be positive")
  }
  
  # Non-linear reduction
  if (!is.null(nonlin.reduction) && !nonlin.reduction %in% c("umap", "tsne")) {
    stop("Invalid 'nonlin.reduction'. Must be one of: 'umap', 'tsne', or NULL")
  }
  if (!is.null(nonlin.reduction) && (!is.numeric(nonlin.dims) || nonlin.dims < 1 || nonlin.dims %% 1 != 0)) {
    stop("'nonlin.dims' must be a positive integer")
  }
  
  # --------------------------
  # Validate input object
  # --------------------------
  if (!inherits(srt, "Seurat")) stop("'srt' must be a Seurat object")
  if (!batch %in% colnames(srt@meta.data)) stop("Batch column '", batch, "' not found in metadata")
  
  assay <- assay %||% DefaultAssay(srt)
  if (!assay %in% Assays(srt)) stop("Assay '", assay, "' not found")
  
  # --------------------------
  # Batch information
  # --------------------------
  batches <- unique(srt@meta.data[[batch]])
  batch.sizes <- table(srt@meta.data[[batch]])
  
  if (verbose) {
    message("Detected ", length(batches), " batches with sizes: ",
            paste(names(batch.sizes), batch.sizes, sep = "=", collapse = ", "))
  }
  
  # Check for very small batches - now stop if too small
  small.batches <- names(batch.sizes)[batch.sizes < bbknn.n_neighbors_within_batch]
  if (length(small.batches) > 0) {
    stop("Batches smaller than n_neighbors_within_batch (", bbknn.n_neighbors_within_batch, "): ",
         paste(small.batches, collapse = ", "), ". Increase batch size or reduce n_neighbors_within_batch to proceed.")
  }
  
  # --------------------------
  # Batch-specific preprocessing
  # --------------------------
  if (verbose) message("\n[1/6] Performing batch-specific preprocessing")
  processed.batches <- list()
  
  for (b in batches) {
    if (verbose) message("  Processing batch: ", b)
    
    # Subset to current batch
    srt.batch <- subset(srt, cells = which(srt@meta.data[[batch]] == b))
    DefaultAssay(srt.batch) <- assay
    
    # Normalization with parameters
    if (recompute.normalization || is.null(srt.batch@assays[[assay]]@data)) {
      if (verbose) message("    Running normalization (", normalization.method, ")")
      srt.batch <- NormalizeData(
        object = srt.batch,
        assay = assay,
        normalization.method = normalization.method,
        scale.factor = scale.factor,
        verbose = FALSE
      )
    } else if (verbose) {
      message("    Using existing normalized data")
    }
    
    # Variable feature selection with parameters
    if (recompute.var.features || length(VariableFeatures(srt.batch, assay = assay)) == 0) {
      if (verbose) message("    Selecting ", var.features.n, " variable features (", var.features.method, ")")
      srt.batch <- FindVariableFeatures(
        object = srt.batch,
        assay = assay,
        selection.method = var.features.method,
        nfeatures = var.features.n,
        verbose = FALSE
      )
      if (verbose) message("    Found ", length(VariableFeatures(srt.batch, assay = assay)), " variable features")
    } else if (verbose) {
      message("    Using existing variable features (", length(VariableFeatures(srt.batch, assay = assay)), ")")
    }
    
    processed.batches[[b]] <- srt.batch
  }
  
  # --------------------------
  # Merge batches and scale
  # --------------------------
  if (verbose) message("\n[2/6] Merging processed batches")
  srt.merged <- merge(
    x = processed.batches[[1]],
    y = processed.batches[-1],
    merge.data = TRUE
  )
  DefaultAssay(srt.merged) <- assay
  
  # Combine variable features (union)
  all.hvgs <- unique(unlist(lapply(processed.batches, function(x) {
    VariableFeatures(x, assay = assay)
  })))
  VariableFeatures(srt.merged, assay = assay) <- all.hvgs
  if (verbose) message("  Combined variable features from all batches: ", length(all.hvgs))
  
  # Scale data
  if (verbose) message("  Scaling data using combined variable features")
  srt.merged <- ScaleData(
    object = srt.merged,
    assay = assay,
    features = all.hvgs,
    verbose = FALSE
  )
  
  # --------------------------
  # Linear dimensionality reduction
  # --------------------------
  if (verbose) message("\n[3/6] Performing linear dimensionality reduction (", linear.reduction, ")")
  valid.reductions <- c("pca", "svd", "ica", Reductions(srt.merged))
  if (!linear.reduction %in% valid.reductions) {
    stop("Invalid 'linear.reduction' method. Must be one of: ", paste(valid.reductions, collapse = ", "))
  }
  
  reduction.name <- paste0("BBKNN_", linear.reduction)
  if (!linear.reduction %in% Reductions(srt.merged)) {
    if (verbose) message("  Computing ", linear.reduction, " with ", linear.reduction.dims, " dimensions")
    
    reduction.args <- list(
      object = srt.merged,
      assay = assay,
      features = all.hvgs,
      reduction.name = reduction.name,
      verbose = FALSE,
      seed.use = seed
    )
    
    srt.merged <- switch(linear.reduction,
      pca = do.call(RunPCA, c(reduction.args, list(npcs = linear.reduction.dims))),
      svd = do.call(RunSVD, c(reduction.args, list(n = linear.reduction.dims))),
      ica = do.call(RunICA, c(reduction.args, list(nics = linear.reduction.dims))),
      stop("Unsupported linear reduction: ", linear.reduction)
    )
  } else {
    reduction.name <- linear.reduction
    if (verbose) message("  Using existing ", linear.reduction, " reduction")
  }
  
  # --------------------------
  # Determine dimensions to use
  # --------------------------
  if (is.null(linear.reduction.dims.use)) {
    if (verbose) message("  Estimating optimal number of dimensions")
    
    dim.est <- tryCatch({
      emb <- Embeddings(srt.merged, reduction = reduction.name)
      min(intrinsicDimension::maxLikGlobalDimEst(emb, k = dim.est.k)[["dim.est"]], linear.reduction.dims)
    }, error = function(e) {
      if (verbose) message("  Dimension estimation failed: ", e$message, ". Falling back to elbow method.")
      if (linear.reduction == "pca" && "pca" %in% Reductions(srt.merged)) {
        variance <- srt.merged@reductions$pca@stdev^2
        elbow <- which.max(diff(cumsum(variance)/sum(variance)))
        min(elbow, linear.reduction.dims)
      } else {
        min(30, linear.reduction.dims)
      }
    })
    
    linear.reduction.dims.use <- seq_len(max(10, ceiling(dim.est)))
    if (verbose) message("  Selected ", length(linear.reduction.dims.use), " dimensions for BBKNN")
  } else if (verbose) {
    message("  Using specified dimensions: ", paste(linear.reduction.dims.use, collapse = ", "))
  }
  
  # --------------------------
  # BBKNN integration
  # --------------------------
  if (verbose) message("\n[4/6] Running BBKNN integration")
  
  # Extract embeddings
  emb <- Embeddings(srt.merged, reduction = reduction.name)[, linear.reduction.dims.use, drop = FALSE]
  
  # Prepare parameters
  params <- list(
    data = emb,
    batch_list = srt.merged[[batch, drop = TRUE]],
    n_neighbors = bbknn.n_neighbors,
    n_neighbors_within_batch = bbknn.n_neighbors_within_batch
  )
  params <- modifyList(params, bbknn.params)
  
  # Run BBKNN
  bem <- tryCatch({
    exec(bbknn$bbknn, !!!params)
  }, error = function(e) {
    stop("BBKNN integration failed: ", e$message, 
         "\nCheck: 1) BBKNN version (pip install --upgrade bbknn), 2) Parameter names (use snake_case like 'n_neighbors'), ",
         "3) Batch sizes >= n_neighbors_within_batch, 4) bbknn.params for invalid entries.")
  })
  
  n.neighbors <- bem$parameters$n_neighbors
  if (verbose) message("  BBKNN completed with ", n.neighbors, " neighbors per cell")
  
  # Store BBKNN graph (using Matrix::sparseMatrix for clarity)
  connectivities <- bem$connectivities
  bbknn.graph <- Matrix::sparseMatrix(
    i = connectivities$nonzero()[0] + 1,
    j = connectivities$nonzero()[1] + 1,
    x = connectivities$data,
    dims = connectivities$shape
  )
  rownames(bbknn.graph) <- colnames(bbknn.graph) <- rownames(emb)
  srt.merged@graphs[["BBKNN"]] <- as.Graph(bbknn.graph)
  srt.merged@graphs[["BBKNN"]]@assay.used <- assay
  
  # Store distance matrix
  distances <- bem$distances
  bbknn.dist <- Matrix::sparseMatrix(
    i = distances$nonzero()[0] + 1,
    j = distances$nonzero()[1] + 1,
    x = distances$data,
    dims = distances$shape
  )
  rownames(bbknn.dist) <- colnames(bbknn.dist) <- rownames(emb)
  srt.merged@graphs[["BBKNN_dist"]] <- as.Graph(bbknn.dist)
  srt.merged@graphs[["BBKNN_dist"]]@assay.used <- assay
  
  # Prepare neighbor object
  val <- split(bbknn.dist@x, rep(1:ncol(bbknn.dist), diff(bbknn.dist@p)))
  pos <- split(bbknn.dist@i + 1, rep(1:ncol(bbknn.dist), diff(bbknn.dist@p)))
  
  has_padding <- FALSE
  idx <- t(mapply(function(x, y) {
    ord <- order(x, decreasing = FALSE)
    out <- y[head(ord, n.neighbors - 1)]
    if (length(out) < n.neighbors - 1) {
      has_padding <<- TRUE
      out <- c(out, rep(y[ord[1]], n.neighbors - 1 - length(out)))
    }
    out
  }, x = val, y = pos))
  idx <- cbind(seq_len(nrow(idx)), idx)
  
  dist <- t(mapply(function(x) {
    ord <- order(x, decreasing = FALSE)
    out <- x[head(ord, n.neighbors - 1)]
    if (length(out) < n.neighbors - 1) {
      out <- c(out, rep(x[ord[1]], n.neighbors - 1 - length(out)))
    }
    c(0, out)
  }, x = val))
  
  if (has_padding && verbose) {
    warning("Some cells have fewer than ", n.neighbors - 1, " neighbors; padded with duplicates. Consider adjusting parameters.")
  }
  
  srt.merged[["BBKNN_neighbors"]] <- new(
    Class = "Neighbor",
    nn.idx = idx,
    nn.dist = dist,
    alg.info = list(),
    cell.names = rownames(emb)
  )
  
  # --------------------------
  # Clustering
  # --------------------------
  if (verbose) message("\n[5/6] Performing clustering with ", cluster.algorithm)
  cluster.algorithm.index <- switch(tolower(cluster.algorithm),
                                    louvain = 1,
                                    slm = 3,
                                    leiden = 4,
                                    stop("Invalid clustering algorithm: ", cluster.algorithm)
  )
  
  srt.merged <- tryCatch({
    srt.merged <- FindClusters(
      object = srt.merged,
      graph.name = "BBKNN",
      resolution = cluster.resolution,
      algorithm = cluster.algorithm.index,
      method = "igraph",
      verbose = FALSE
    )
    srt.merged[["BBKNNclusters"]] <- Idents(srt.merged)
    if (verbose) message("  Clustering completed with ", length(unique(Idents(srt.merged))), " clusters")
    srt.merged
  }, error = function(e) {
    warning("Clustering failed: ", e$message, ". Results not stored.")
    srt.merged
  })
  
  # --------------------------
  # Non-linear dimensionality reduction
  # --------------------------
  if (!is.null(nonlin.reduction)) {
    if (verbose) message("\n[6/6] Running ", nonlin.reduction, " for visualization")
    reduction.name.nonlin <- paste0("BBKNN_", nonlin.reduction)
    
    reduction.args <- list(
      object = srt.merged,
      reduction = reduction.name,
      dims = linear.reduction.dims.use,
      seed.use = seed,
      verbose = FALSE
    )
    
    if (nonlin.reduction == "umap") {
      reduction.args <- c(reduction.args, list(
        min.dist = umap.min.dist,
        n.neighbors = umap.n.neighbors
      ))
    }
    
    srt.merged <- if (nonlin.reduction == "umap") {
      do.call(RunUMAP, c(reduction.args, list(reduction.name = reduction.name.nonlin)))
    } else {
      do.call(RunTSNE, c(reduction.args, list(reduction.name = reduction.name.nonlin)))
    }
    
    if (verbose) message("  ", nonlin.reduction, " results stored as '", reduction.name.nonlin, "'")
  }
  
  # --------------------------
  # Final setup
  # --------------------------
  DefaultAssay(srt.merged) <- assay
  if (verbose) message("\nBBKNN pipeline completed successfully")
  
  return(srt.merged)
}