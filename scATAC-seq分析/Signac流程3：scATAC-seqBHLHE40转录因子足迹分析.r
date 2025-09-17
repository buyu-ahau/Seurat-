# BHLHE40转录因子足迹分析
# 当前日期：2025-09-03
# 用户：buyu-ahau

# 加载必要的包
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Rsamtools)

# 设置输出目录
main_output_dir <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/footprinting_analysis_results"
bhlhe40_output_dir <- file.path(main_output_dir, "BHLHE40足迹分析")
dir.create(bhlhe40_output_dir, recursive = TRUE, showWarnings = FALSE)

# 创建日志文件
log_file <- file.path(bhlhe40_output_dir, "BHLHE40_footprint_analysis_log.txt")
sink(log_file)
cat("BHLHE40转录因子足迹分析\n")
cat("=======================\n")
cat("日期: 2025-09-03\n")
cat("用户: buyu-ahau\n\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# 加载Seurat对象
cat("加载Seurat对象...\n")
atac_obj_file <- file.path(main_output_dir, "atac_with_updated_fragments.rds")
atac_obj <- readRDS(atac_obj_file)
cat("Seurat对象加载完成\n\n")

# 获取可用的细胞类型
all_cell_types <- unique(atac_obj$seurat_clusters)
cat("可用的细胞类型:\n")
print(all_cell_types)
cat("\n")

# 检查Seurat对象中是否包含BHLHE40基序
cat("检查BHLHE40基序信息...\n")
motif_obj <- Motifs(atac_obj[["ATAC"]])
if(!is.null(motif_obj) && !is.null(motif_obj@pwm)) {
  cat("对象包含", length(motif_obj@pwm), "个基序\n")
  
  # 查找BHLHE40基序
  bhlhe40_motif <- "MA0464.2"  # BHLHE40的基序ID
  if(bhlhe40_motif %in% names(motif_obj@pwm)) {
    cat("找到BHLHE40基序ID:", bhlhe40_motif, "\n")
  } else {
    cat("警告: 未找到BHLHE40基序ID:", bhlhe40_motif, "\n")
    cat("尝试搜索包含'BHLHE40'的基序...\n")
    
    bhlhe40_candidates <- grep("BHLHE40", names(motif_obj@pwm), value = TRUE)
    if(length(bhlhe40_candidates) > 0) {
      cat("找到以下可能的BHLHE40基序:", paste(bhlhe40_candidates, collapse = ", "), "\n")
      bhlhe40_motif <- bhlhe40_candidates[1]
      cat("将使用:", bhlhe40_motif, "\n")
    } else {
      cat("错误: 未找到BHLHE40相关基序\n")
      stop("未找到BHLHE40基序")
    }
  }
} else {
  cat("错误: Seurat对象中没有motif信息\n")
  stop("Seurat对象中没有motif信息")
}
cat("\n")

# 获取猪基因组文件
cat("设置猪基因组参考...\n")
genome_path <- "/disk192/users_dir/buyu/2.参考基因组/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
if(file.exists(genome_path)) {
  pig_genome <- Rsamtools::FaFile(genome_path)
  cat("猪基因组文件加载成功\n")
} else {
  cat("错误: 猪基因组文件不存在:", genome_path, "\n")
  stop("猪基因组文件不存在")
}
cat("\n")

# 执行BHLHE40足迹分析
cat("开始BHLHE40足迹分析...\n")
cat("基序ID:", bhlhe40_motif, "\n")
start_time <- Sys.time()

# 对所有细胞类型进行足迹分析
cat("\n分析所有细胞类型...\n")

# 执行足迹分析
cat("执行足迹分析...\n")
tryCatch({
  footprint_data <- Footprint(
    object = atac_obj,
    motif.name = bhlhe40_motif,
    assay = "ATAC", 
    group.by = "seurat_clusters",
    idents = all_cell_types,  # 使用所有细胞类型
    upstream = 500,
    downstream = 500,
    genome = pig_genome
  )
  
  cat("足迹分析成功!\n")
  
  # 保存足迹数据
  result_file <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_all_cells.rds")
  saveRDS(footprint_data, result_file)
  cat("结果保存至:", result_file, "\n")
  
  # 可视化足迹数据 - 所有细胞类型在同一图中
  cat("生成所有细胞类型在同一图中的可视化...\n")
  
  # 单一图中显示所有细胞类型
  p_all_in_one <- PlotFootprint(
    object = footprint_data,
    features = bhlhe40_motif,
    idents = all_cell_types,  # 包含所有细胞类型
    group.by = "seurat_clusters"  # 按细胞类型分组
  ) + 
  labs(title = "BHLHE40 Footprint - All Cell Types") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
  
  # 保存所有细胞类型在同一图中的足迹图
  pdf_file_all <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_all_cells_in_one_plot.pdf")
  pdf(pdf_file_all, width = 12, height = 10)
  print(p_all_in_one)
  dev.off()
  
  png_file_all <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_all_cells_in_one_plot.png")
  png(png_file_all, width = 1800, height = 1500, res = 150)
  print(p_all_in_one)
  dev.off()
  
  cat("所有细胞类型在同一图中的足迹图已保存至:", pdf_file_all, "和", png_file_all, "\n")
  
  # 创建免疫细胞和非免疫细胞的分组图
  cat("创建免疫细胞和非免疫细胞的分组足迹图...\n")
  
  # 免疫细胞组
  immune_cells <- c("B_cells", "CD4_CD8_T_cells", "Macrophages", "NK_cells", "Neutrophils", "Plasma_cells")
  # 非免疫细胞组
  non_immune_cells <- c("Endothelial_cells", "Epithelial_cells", "Fibroblasts")
  
  # 为免疫细胞绘制单一图
  p_immune <- PlotFootprint(
    object = footprint_data,
    features = bhlhe40_motif,
    idents = immune_cells,
    group.by = "seurat_clusters"
  ) + 
  labs(title = "BHLHE40 Footprint - Immune Cells") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
  
  # 为非免疫细胞绘制单一图
  p_non_immune <- PlotFootprint(
    object = footprint_data,
    features = bhlhe40_motif,
    idents = non_immune_cells,
    group.by = "seurat_clusters"
  ) + 
  labs(title = "BHLHE40 Footprint - Non-Immune Cells") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
  
  # 保存免疫细胞和非免疫细胞的足迹图
  pdf_file_immune <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_immune_cells.pdf")
  pdf_file_non_immune <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_non_immune_cells.pdf")
  
  pdf(pdf_file_immune, width = 12, height = 10)
  print(p_immune)
  dev.off()
  
  pdf(pdf_file_non_immune, width = 12, height = 10)
  print(p_non_immune)
  dev.off()
  
  png_file_immune <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_immune_cells.png")
  png_file_non_immune <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_non_immune_cells.png")
  
  png(png_file_immune, width = 1800, height = 1500, res = 150)
  print(p_immune)
  dev.off()
  
  png(png_file_non_immune, width = 1800, height = 1500, res = 150)
  print(p_non_immune)
  dev.off()
  
  cat("免疫细胞和非免疫细胞的足迹图已保存\n")
  
  # 尝试创建一个自定义版本的足迹图，优化显示效果
  cat("创建优化版本的足迹图...\n")
  
  # 尝试使用更多颜色和更清晰的线条
  custom_colors <- c(
    "B_cells" = "#E41A1C",            # 红色
    "CD4_CD8_T_cells" = "#377EB8",    # 蓝色
    "Plasma_cells" = "#4DAF4A",       # 绿色
    "NK_cells" = "#984EA3",           # 紫色
    "Endothelial_cells" = "#FF7F00",  # 橙色
    "Neutrophils" = "#FFFF33",        # 黄色
    "Macrophages" = "#A65628",        # 棕色
    "Epithelial_cells" = "#F781BF",   # 粉色
    "Fibroblasts" = "#999999"         # 灰色
  )
  
  # 创建自定义足迹图
  p_custom <- PlotFootprint(
    object = footprint_data,
    features = bhlhe40_motif,
    idents = all_cell_types,
    group.by = "seurat_clusters"
  ) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "BHLHE40 Footprint Analysis Across All Cell Types",
    subtitle = paste("Motif:", bhlhe40_motif, "- Sus scrofa (Pig) Genome"),
    x = "Distance from motif (bp)",
    y = "Tn5 insertion profile"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
  
  # 保存自定义足迹图
  pdf_file_custom <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_all_cells_custom.pdf")
  pdf(pdf_file_custom, width = 14, height = 10)
  print(p_custom)
  dev.off()
  
  png_file_custom <- file.path(bhlhe40_output_dir, "footprint_BHLHE40_all_cells_custom.png")
  png(png_file_custom, width = 2000, height = 1500, res = 150)
  print(p_custom)
  dev.off()
  
  cat("优化版本的足迹图已保存\n")
  
}, error = function(e) {
  cat("错误: 足迹分析失败:", e$message, "\n")
})
# 分析完成
end_time <- Sys.time()
cat("\n\nBHLHE40足迹分析完成\n")
cat("总耗时:", difftime(end_time, start_time, units = "mins"), "分钟\n")
# 保存当前会话信息
cat("\n保存会话信息...\n")
session_info_file <- file.path(bhlhe40_output_dir, "session_info_BHLHE40.txt")
writeLines(capture.output(sessionInfo()), session_info_file)
cat("会话信息已保存至:", session_info_file, "\n")
cat("\n=== 分析完成 ===\n")
cat("结果保存至:", bhlhe40_output_dir, "\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()
# 打印分析完成信息
cat("\n=== BHLHE40足迹分析完成 ===\n")
cat("结果保存至:", bhlhe40_output_dir, "\n")
cat("日志文件:", log_file, "\n")





