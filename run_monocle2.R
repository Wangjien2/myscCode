#!/usr/bin/env Rscript

## 加载必要的库
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(DDRTree)
  library(BiocGenerics)
  library(Biobase)
  library(VGAM)
  library(optparse)
  library(rlang)
})

## 定义辅助函数：检查R包是否安装
check_R <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("缺少必要的R包: ", paste(missing, collapse = ", "))
  }
}

## 1. 从Seurat转为Monocle2 CellDataSet
prepare_cds_from_seurat <- function(srt, assay = NULL, slot = "counts", expressionFamily = "negbinomial.size") {
  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- as.sparse(GetAssayData(srt, assay = assay, slot = slot))
  p_data <- srt@meta.data
  f_data <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))
  
  cds <- monocle::newCellDataSet(
    expr_matrix,
    phenoData   = new("AnnotatedDataFrame", data = p_data),
    featureData = new("AnnotatedDataFrame", data = f_data),
    expressionFamily = do.call(get(expressionFamily, envir = getNamespace("VGAM")), args = list())
  )
  
  if (any(c("negbinomial", "negbinomial.size") %in% expressionFamily)) {
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- suppressWarnings(suppressMessages(BiocGenerics::estimateDispersions(cds)))
  }
  
  return(cds)
}

## 2. 特征选择
select_features <- function(srt, cds, assay = NULL, feature_type = "HVF", disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit") {
  features <- NULL
  assay <- assay %||% DefaultAssay(srt)
  
  if (feature_type == "HVF") {
    features <- VariableFeatures(srt, assay = assay)
    if (length(features) == 0) {
      srt <- FindVariableFeatures(srt, assay = assay)
      features <- VariableFeatures(srt, assay = assay)
    }
  } else if (feature_type == "Disp") {
    features <- subset(monocle::dispersionTable(cds), eval(rlang::parse_expr(disp_filter)))$gene_id
  }
  
  message("特征数量: ", length(features))
  cds <- monocle::setOrderingFilter(cds, features)
  return(list(cds = cds, features = features, srt = srt))
}

## 3. 运行Monocle2降维和轨迹推断
build_monocle2 <- function(cds, max_components = 2, reduction_method = "DDRTree", 
                          norm_method = "log", residualModelFormulaStr = NULL, 
                          pseudo_expr = 1, root_state = NULL) {
  cds <- monocle::reduceDimension(
    cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr
  )
  cds <- monocle::orderCells(cds, root_state = root_state)
  return(cds)
}

## 4. 绘制轨迹
plot_monocle2_trajectory <- function(srt, cds, reduction_method = "DDRTree", 
                                    annotation = NULL, output_dir) {
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 绘制状态轨迹图
  p1 <- DimPlot(srt, group.by = "Monocle2_State", reduction = reduction_method, 
               label = TRUE, force = TRUE) + 
    ggtitle("Monocle2 细胞状态轨迹")
  ggsave(file.path(output_dir, "monocle2_state_trajectory.png"), 
         plot = p1, width = 10, height = 8, dpi = 300)
  
  # 绘制拟时序轨迹图
  p2 <- FeaturePlot(srt, features = "Monocle2_Pseudotime", reduction = reduction_method) +
    ggtitle("Monocle2 拟时序轨迹")
  ggsave(file.path(output_dir, "monocle2_pseudotime_trajectory.png"), 
         plot = p2, width = 10, height = 8, dpi = 300)
  
  # 绘制注释轨迹图（如果提供）
  if (!is.null(annotation)) {
    p3 <- DimPlot(srt, group.by = annotation, reduction = reduction_method, 
                 label = TRUE, force = TRUE) +
      ggtitle(paste("按", annotation, "分组的轨迹"))
    ggsave(file.path(output_dir, paste0("monocle2_", annotation, "_trajectory.png")), 
           plot = p3, width = 10, height = 8, dpi = 300)
  }
}

## 5. 主函数
RunMonocle2Trajectory <- function(
    srt, 
    annotation = NULL, 
    assay = NULL, 
    slot = "counts", 
    expressionFamily = "negbinomial.size", 
    features = NULL, 
    feature_type = "HVF", 
    disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit", 
    max_components = 2, 
    reduction_method = "DDRTree", 
    norm_method = "log", 
    residualModelFormulaStr = NULL, 
    pseudo_expr = 1, 
    root_state = NULL, 
    seed = 11,
    plot = TRUE,
    output_dir = "./monocle2_plots"
) {
  set.seed(seed)
  check_R(c("monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM", "Seurat", "ggplot2"))
  
  ## Step 1: 转换 Seurat -> CDS
  cds <- prepare_cds_from_seurat(srt, assay, slot, expressionFamily)
  
  ## Step 2: 特征选择
  if (is.null(features)) {
    res <- select_features(srt, cds, assay, feature_type, disp_filter)
    cds <- res$cds
    features <- res$features
    srt <- res$srt
  } else {
    cds <- monocle::setOrderingFilter(cds, features)
  }
  
  ## Step 3: 构建 Monocle2
  cds <- build_monocle2(cds, max_components, reduction_method, norm_method, 
                       residualModelFormulaStr, pseudo_expr, root_state)
  
  ## Step 4: 保存结果到 Seurat
  embeddings <- t(cds@reducedDimS)
  colnames(embeddings) <- paste0(cds@dim_reduce_type, "_", 1:ncol(embeddings))
  srt[[cds@dim_reduce_type]] <- CreateDimReducObject(
    embeddings = embeddings, 
    key = paste0(cds@dim_reduce_type, "_"), 
    assay = assay
  )
  srt[["Monocle2_State"]] <- cds[["State"]]
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(cds = cds, features = features)
  
  ## Step 5: 绘图
  if (plot) {
    plot_monocle2_trajectory(srt, cds, reduction_method, annotation, output_dir)
  }
  
  return(srt)
}

## 解析命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              help = "输入的Seurat对象RDS文件路径", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", 
              help = "输出的Seurat对象RDS文件路径", metavar = "FILE"),
  make_option(c("-p", "--plot_dir"), type = "character", default = "./monocle2_plots",
              help = "轨迹图输出目录 [default: %default]"),
  make_option(c("-a", "--assay"), type = "character", default = NULL,
              help = "使用的Seurat assay名称"),
  make_option(c("-s", "--slot"), type = "character", default = "counts",
              help = "使用的数据槽 [default: %default]"),
  make_option(c("-f", "--feature_type"), type = "character", default = "HVF",
              help = "特征选择类型 (HVF/Disp) [default: %default]"),
  make_option(c("-r", "--reduction"), type = "character", default = "DDRTree",
              help = "降维方法 [default: %default]"),
  make_option(c("-c", "--components"), type = "integer", default = 2,
              help = "降维维度 [default: %default]"),
  make_option(c("-t", "--annotation"), type = "character", default = NULL,
              help = "用于注释的元数据列名"),
  make_option(c("-root", "--root_state"), type = "integer", default = NULL,
              help = "根状态ID"),
  make_option(c("-seed", "--seed"), type = "integer", default = 11,
              help = "随机种子 [default: %default]"),
  make_option(c("-noplot", "--no_plot"), action = "store_false", default = TRUE,
              help = "不生成轨迹图")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

## 参数验证
if (is.null(args$input) || is.null(args$output)) {
  print_help(parser)
  stop("必须提供输入和输出文件路径", call. = FALSE)
}

## 主流程执行
message("加载Seurat对象: ", args$input)
srt <- readRDS(args$input)

message("运行Monocle2轨迹分析...")
srt <- RunMonocle2Trajectory(
  srt = srt,
  annotation = args$annotation,
  assay = args$assay,
  slot = args$slot,
  feature_type = args$feature_type,
  max_components = args$components,
  reduction_method = args$reduction,
  root_state = args$root_state,
  seed = args$seed,
  plot = args$no_plot,
  output_dir = args$plot_dir
)

message("保存结果到: ", args$output)
saveRDS(srt, file = args$output)

message("分析完成!")
