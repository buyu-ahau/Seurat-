library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Rsamtools)
# 设置目录
main_dir <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/footprinting_analysis_results"
output_dir <- file.path(main_dir, "BHLHE40足迹分析")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# 加载数据
atac_obj <- readRDS(file.path(main_dir, "atac_with_updated_fragments.rds"))
all_cell_types <- unique(atac_obj$seurat_clusters)
# 获取猪基因组
genome_path <- "/disk192/users_dir/buyu/2.参考基因组/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
pig_genome <- Rsamtools::FaFile(genome_path)
# 定义转录因子列表及其JASPAR ID
tf_list <- c(
  "PAX5" = "MA0014.3",   # PAX5
  "GATA3" = "MA0037.3",  # GATA3
  "SPI1" = "MA0080.5",   # SPI1 (PU.1)
  "EOMES" = "MA0800.1",  # EOMES
  "CEBPE" = "MA0837.1",  # CEBPE
  "IRF4" = "MA0050.2"    # IRF4
)
# 检查转录因子是否存在
motif_obj <- Motifs(atac_obj[["ATAC"]])
available_motifs <- names(motif_obj@pwm)
available_tfs <- intersect(tf_list, available_motifs)
if(length(available_tfs) == 0) {
  stop("没有找到任何指定的转录因子")
} else {
  print(paste("找到", length(available_tfs), "个转录因子:", paste(names(available_tfs), collapse=", ")))
}
# 定义细胞组
cell_groups <- list(
  "所有细胞" = all_cell_types,
  "免疫细胞" = c("B_cells", "CD4_CD8_T_cells", "Macrophages", "NK_cells", "Neutrophils", "Plasma_cells"),
  "非免疫细胞" = c("Endothelial_cells", "Epithelial_cells", "Fibroblasts")
)
# 执行足迹分析
all_results <- list()
for(tf_name in names(tf_list)) {
  tf_id <- tf_list[tf_name]
  
  if(tf_id %in% available_motifs) {
    print(paste("分析", tf_name, "(", tf_id, ")"))
    
    # 创建TF特定目录
    tf_dir <- file.path(output_dir, tf_name)
    dir.create(tf_dir, showWarnings = FALSE)
    
    # 对所有细胞分析
    footprint_data <- Footprint(
      object = atac_obj,
      motif.name = tf_id,
      assay = "ATAC",
      group.by = "seurat_clusters",
      idents = all_cell_types,
      upstream = 500,
      downstream = 500,
      genome = pig_genome
    )
    
    # 保存结果
    saveRDS(footprint_data, file.path(tf_dir, paste0(tf_name, "_footprint_all_cells.rds")))
    all_results[[tf_name]] <- footprint_data
    
    # 创建可视化
    # 所有细胞在一张图
    p_all <- PlotFootprint(
      object = footprint_data,
      features = tf_id,
      idents = all_cell_types,
      group.by = "seurat_clusters"
    ) + 
    labs(title = paste(tf_name, "Footprint - All Cell Types"))
    
    pdf(file.path(tf_dir, paste0(tf_name, "_all_cells.pdf")), width = 12, height = 10)
    print(p_all)
    dev.off()
    
    # 免疫细胞图
    p_immune <- PlotFootprint(
      object = footprint_data,
      features = tf_id,
      idents = cell_groups[["免疫细胞"]],
      group.by = "seurat_clusters"
    ) + 
    labs(title = paste(tf_name, "Footprint - Immune Cells"))
    
    pdf(file.path(tf_dir, paste0(tf_name, "_immune_cells.pdf")), width = 12, height = 8)
    print(p_immune)
    dev.off()
    
    # 非免疫细胞图
    p_non_immune <- PlotFootprint(
      object = footprint_data,
      features = tf_id,
      idents = cell_groups[["非免疫细胞"]],
      group.by = "seurat_clusters"
    ) + 
    labs(title = paste(tf_name, "Footprint - Non-Immune Cells"))
    
    pdf(file.path(tf_dir, paste0(tf_name, "_non_immune_cells.pdf")), width = 12, height = 8)
    print(p_non_immune)
    dev.off()
  } else {
    print(paste("警告: 未找到", tf_name, "(", tf_id, ")"))
  }
}





###scRNA-seq数据转录因子分析结果footprint分析
# 转录因子足迹分析 - 使用猪基因组
# 作者: buyu-ahau
# 日期: 2025-09-03

# 加载必要的包
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Rsamtools)  # 用于加载基因组
library(BSgenome)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(pheatmap)
library(RColorBrewer)

# 设置输出目录
main_dir <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/footprinting_analysis_results"
output_dir <- file.path(main_dir, "TF_footprints_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 创建日志文件
log_file <- file.path(output_dir, "footprint_analysis_log.txt")
sink(log_file, split = TRUE) # split = TRUE 同时输出到控制台和日志文件

cat("转录因子足迹分析 - 使用猪基因组\n")
cat("====================================\n")
cat("日期: 2025-09-03\n")
cat("用户: buyu-ahau\n\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# 加载数据
cat("加载Seurat对象...\n")
atac_obj_file <- file.path(main_dir, "atac_with_updated_fragments.rds")
atac_obj <- readRDS(atac_obj_file)
cat("Seurat对象加载完成\n\n")

# 加载猪基因组
cat("加载猪基因组...\n")
genome_path <- "/disk192/users_dir/buyu/2.参考基因组/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
pig_genome <- FaFile(genome_path)
cat("猪基因组加载完成\n\n")

# 获取可用的细胞类型
all_cell_types <- unique(atac_obj$seurat_clusters)
cat("可用的细胞类型:\n")
print(all_cell_types)
cat("\n")

# 检查Seurat对象中已有的motifs
cat("检查Seurat对象中已有的motifs...\n")
motif_obj <- Motifs(atac_obj[["ATAC"]])
if(!is.null(motif_obj) && !is.null(motif_obj@pwm)) {
  cat("对象包含", length(motif_obj@pwm), "个基序\n")
} else {
  cat("警告: Seurat对象中没有motif信息，将尝试获取JASPAR基序\n")
}
cat("\n")

# 定义31个转录因子列表
tf_list <- c(
  "FOXP1" = "MA0481.3",
  "TBX21" = "MA0809.2",
  "IKZF2" = "MA0486.2",
  "IRF2" = "MA0051.1",
  "ZEB1" = "MA0103.3",
  "GFI1" = "MA0038.2",
  "RORA" = "MA0071.1",
  "TCF7" = "MA0769.2",
  "KLF2" = "MA0599.1",
  "EOMES" = "MA0800.1",
  "E2F1" = "MA0024.3",
  "STAT1" = "MA0137.3",
  "FOXO1" = "MA0480.1",
  "RUNX3" = "MA0687.1",
  "GATA3" = "MA0037.3",
  "ETS1" = "MA0098.3",
  "FOXP3" = "MA0481.3",
  "FLI1" = "MA0475.2",
  "IRF4" = "MA0050.2",
  "PRDM1" = "MA0508.3",
  "XBP1" = "MA0839.1",
  "EBF1" = "MA0154.4",
  "PAX5" = "MA0014.3",
  "PPARG" = "MA0066.1",
  "SPI1" = "MA0080.5",
  "STAT3" = "MA0144.2",
  "NFKB1" = "MA0105.4",
  "RELA" = "MA0107.1",
  "IRF8" = "MA0652.1",
  "CEBPB" = "MA0466.2",
  "SOX9" = "MA0077.1",
  "BHLHE40" = "MA0464.2"
)

# 检查哪些转录因子ID存在
available_motifs <- names(motif_obj@pwm)
available_tfs <- tf_list[tf_list %in% available_motifs]

cat("找到", length(available_tfs), "/", length(tf_list), "个有效的转录因子ID\n")
cat("将分析以下转录因子:\n")
for(tf_name in names(available_tfs)) {
  cat("  ", tf_name, ":", available_tfs[tf_name], "\n")
}
cat("\n")

# 定义细胞组
cell_groups <- list(
  "所有细胞" = all_cell_types,
  "免疫细胞" = c("B_cells", "CD4_CD8_T_cells", "Macrophages", "NK_cells", "Neutrophils", "Plasma_cells"),
  "非免疫细胞" = c("Endothelial_cells", "Epithelial_cells", "Fibroblasts")
)

# 定义自定义颜色
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

# 设置进度跟踪
start_time <- Sys.time()
total_tfs <- length(names(tf_list))
current_tf <- 0

# 执行足迹分析
cat("\n开始足迹分析...\n")
all_results <- list()

for(tf_name in names(tf_list)) {
  current_tf <- current_tf + 1
  tf_id <- tf_list[tf_name]
  
  cat(sprintf("[%d/%d] 分析 %s (%s)\n", current_tf, total_tfs, tf_name, tf_id))
  
  # 只处理可用的转录因子
  if(tf_id %in% available_motifs) {
    # 创建TF特定目录
    tf_dir <- file.path(output_dir, tf_name)
    dir.create(tf_dir, showWarnings = FALSE)
    
    # 尝试执行足迹分析
    tryCatch({
      # 对所有细胞分析
      footprint_data <- Footprint(
        object = atac_obj,
        motif.name = tf_id,
        assay = "ATAC",
        group.by = "seurat_clusters",
        idents = all_cell_types,
        upstream = 500,
        downstream = 500,
        genome = pig_genome
      )
      
      # 保存结果
      result_file <- file.path(tf_dir, paste0(tf_name, "_footprint_all_cells.rds"))
      saveRDS(footprint_data, result_file)
      all_results[[tf_name]] <- footprint_data
      cat("  结果已保存至:", result_file, "\n")
      
      # 创建可视化 - 所有细胞在一张图
      p_all <- PlotFootprint(
        object = footprint_data,
        features = tf_id,
        idents = all_cell_types,
        group.by = "seurat_clusters"
      ) + 
      labs(title = paste(tf_name, "Footprint - All Cell Types"))
      
      pdf_file <- file.path(tf_dir, paste0(tf_name, "_all_cells.pdf"))
      pdf(pdf_file, width = 12, height = 10)
      print(p_all)
      dev.off()
      
      # 同时保存PNG格式以便快速查看
      png_file <- file.path(tf_dir, paste0(tf_name, "_all_cells.png"))
      png(png_file, width = 1800, height = 1500, res = 150)
      print(p_all)
      dev.off()
      
      cat("  所有细胞图表已保存\n")
      
      # 免疫细胞图
      p_immune <- PlotFootprint(
        object = footprint_data,
        features = tf_id,
        idents = cell_groups[["免疫细胞"]],
        group.by = "seurat_clusters"
      ) + 
      labs(title = paste(tf_name, "Footprint - Immune Cells"))
      
      pdf_file_immune <- file.path(tf_dir, paste0(tf_name, "_immune_cells.pdf"))
      pdf(pdf_file_immune, width = 12, height = 8)
      print(p_immune)
      dev.off()
      
      # 非免疫细胞图
      p_non_immune <- PlotFootprint(
        object = footprint_data,
        features = tf_id,
        idents = cell_groups[["非免疫细胞"]],
        group.by = "seurat_clusters"
      ) + 
      labs(title = paste(tf_name, "Footprint - Non-Immune Cells"))
      
      pdf_file_non_immune <- file.path(tf_dir, paste0(tf_name, "_non_immune_cells.pdf"))
      pdf(pdf_file_non_immune, width = 12, height = 8)
      print(p_non_immune)
      dev.off()
      
      cat("  免疫细胞和非免疫细胞图表已保存\n")
      
      # 创建自定义足迹图 - 使用更好的颜色
      p_custom <- PlotFootprint(
        object = footprint_data,
        features = tf_id,
        idents = all_cell_types,
        group.by = "seurat_clusters"
      ) +
      scale_color_manual(values = custom_colors) +
      labs(
        title = paste(tf_name, "Footprint Analysis Across All Cell Types"),
        subtitle = "Sus scrofa (Pig) scATAC-seq data",
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
      
      pdf_file_custom <- file.path(tf_dir, paste0(tf_name, "_all_cells_custom.pdf"))
      pdf(pdf_file_custom, width = 14, height = 10)
      print(p_custom)
      dev.off()
      
      cat("  自定义图表已保存\n")
      
    }, error = function(e) {
      cat("  错误: 足迹分析失败:", conditionMessage(e), "\n")
    })
  } else {
    cat("  警告: 未找到", tf_name, "(", tf_id, ")\n")
  }
}

# 尝试创建热图
if(length(all_results) >= 5) {
  cat("\n创建转录因子富集热图...\n")
  tryCatch({
    # 准备热图数据
    enrichment_matrix <- matrix(0, nrow = length(all_results), ncol = length(all_cell_types))
    rownames(enrichment_matrix) <- names(all_results)
    colnames(enrichment_matrix) <- all_cell_types
    
    for(tf_name in names(all_results)) {
      footprint_data <- all_results[[tf_name]]
      enrichment_data <- footprint_data@enrichment
      
      for(cell_type in names(enrichment_data)) {
        if(cell_type %in% colnames(enrichment_matrix)) {
          enrichment_matrix[tf_name, cell_type] <- enrichment_data[[cell_type]]
        }
      }
    }
    
    # 绘制热图
    pdf_file <- file.path(output_dir, "TF_enrichment_heatmap.pdf")
    pdf(pdf_file, width = 12, height = 10)
    
    pheatmap(
      enrichment_matrix,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "Transcription Factor Footprint Enrichment",
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10
    )
    
    dev.off()
    
    # 同时保存PNG格式
    png_file <- file.path(output_dir, "TF_enrichment_heatmap.png")
    png(png_file, width = 1800, height = 1500, res = 150)
    
    pheatmap(
      enrichment_matrix,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "Transcription Factor Footprint Enrichment",
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10
    )
    
    dev.off()
    
    cat("  热图已保存至:", pdf_file, "和", png_file, "\n")
    
  }, error = function(e) {
    cat("  创建热图失败:", conditionMessage(e), "\n")
  })
}

# 分析完成
end_time <- Sys.time()
cat("\n\n足迹分析完成\n")
cat("总共分析了", length(all_results), "个转录因子\n")
cat("总耗时:", difftime(end_time, start_time, units = "mins"), "分钟\n")
cat("结果保存在:", output_dir, "\n")
cat("结束时间:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

# 保存会话信息
cat("\n保存会话信息...\n")
session_info_file <- file.path(output_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_file)
cat("会话信息已保存至:", session_info_file, "\n")

# 关闭日志文件
sink()

# 打印完成信息
cat("\n=== 转录因子足迹分析完成 ===\n")
cat("结果保存至:", output_dir, "\n")
cat("日志文件:", log_file, "\n")