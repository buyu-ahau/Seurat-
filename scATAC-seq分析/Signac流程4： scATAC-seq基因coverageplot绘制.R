# ---适配GTF注释以匹配Seurat对象并绘图 ---
# 目的:
# 1. 遵循最正确的原则：保持原始Seurat对象的坐标系统不变。
# 2. 将从外部GTF文件加载的基因注释，强制转换为与Seurat对象一致的Ensembl命名风格 ('1', 23456'13'...)。
# 3. 最终成功为目标基因BHLHE40绘制覆盖度图。
# ----------------------------------------------------------------
#!/usr/bin/env Rscript

# -----------------------------------------------------------------
# 脚本目的：绘制 BHLHE40 基因的 scATAC-seq 覆盖度图 (ymax=300)
#
# 请确保已激活 conda 环境, e.g., conda activate signac_motif
# -----------------------------------------------------------------

# ===================================================================
# 1. 加载所有需要的 R 包
# ===================================================================
print("--- 正在加载 R 包 ---")
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))
print("R 包加载完毕。")
print("---------------------------------")


# ===================================================================
# 2. scATAC-seq 分析：加载数据并添加注释
# ===================================================================
print("--- 开始 scATAC-seq 绘图部分 ---")

# --- 步骤 2.1: 加载 scATAC-seq Seurat 对象 ---
atac_file_path <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/obj_relinked_annotated.rds"
print(paste("正在从", atac_file_path, "加载 scATAC 对象..."))
obj_atac <- readRDS(atac_file_path)
print("scATAC 对象加载成功。")

# --- 步骤 2.2: 加载并设置基因组注释 (解决 'Gene not found' 错误) ---
print("正在加载 GTF 基因组注释文件...")
gtf_path <- "/disk192/users_dir/buyu/2.参考基因组/Sus_scrofa.Sscrofa11.1.114.gtf.gz"
pig_annotation_from_gtf <- rtracklayer::import(gtf_path)
print("GTF 加载成功。")

print("正在适配 GTF 染色体命名风格 (e.g., '13' vs 'chr13')...")
# 强制 GTF 风格与 Seurat 对象一致 (使用 '1', '2' 而不是 'chr1', 'chr2')
seqlevelsStyle(pig_annotation_from_gtf) <- "Ensembl"

print("正在将注释添加到 Seurat 对象中...")
# 将修复好的注释重新赋给加载的对象
Annotation(obj_atac) <- pig_annotation_from_gtf
print("注释添加成功！")

# ===================================================================
# 3. 绘制 BHLHE40 基因图 (ymax=300)
# ===================================================================
roi_gene <- "BHLHE40"
print(paste("已设定目标基因为:", roi_gene))
print(paste("正在为", roi_gene, "生成覆盖度图 (ymax=300)..."))

p_bhlhe40 <- CoveragePlot(
  object = obj_atac,
  region = roi_gene,
  extend.upstream = 2000,    # 向上游（左侧）扩展 2000 bp
  extend.downstream = 2000,  # 向下游（右侧）扩展 2000 bp
  annotation = TRUE,         # 显示基因结构注释
  peaks = TRUE,              # 显示已识别的 Peaks
  tile = TRUE,               # 显示 Tile 图，展示片段密度
  ymax = 300                 # 按照你的要求设置为 300
)

# --- 步骤 4: 保存图像 ---
print("正在保存 BHLHE40 基因图...")
ggsave("BHLHE40_CoveragePlot_ymax300.png", plot = p_bhlhe40, width = 8, height = 6, dpi = 300)
ggsave("BHLHE40_CoveragePlot_ymax300.pdf", plot = p_bhlhe40, width = 8, height = 6)

print("---------------------------------")
print("--- 任务完成 ---")
print("请检查你的工作目录下的 'BHLHE40_CoveragePlot_ymax300.png' 文件。")




# --- (前面步骤的代码保持不变) ---
#scRNA-seq合并CD4与CD8细胞绘制小提琴图
# --- 步骤 0 到 3: 加载、合并、定义顺序和颜色 ---
library(Seurat)
library(qs)
library(ggplot2)
library(stringr)

qs_file_path <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/1.scRNA-seq数据/rna_final_annotated_updated.qs"
obj_rna <- qread(qs_file_path)

obj_rna$celltype_plot <- as.character(obj_rna$celltype)
obj_rna$celltype_plot[obj_rna$celltype_plot == "CD4 cells" | obj_rna$celltype_plot == "CD8 cells"] <- "CD4_CD8_T_cells"
obj_rna$celltype_plot <- str_replace_all(obj_rna$celltype_plot, " ", "_")

celltype_order <- c("B_cells", "CD4_CD8_T_cells", "DC_cells", "Endothelial_cells", "Epithelial_cells", "Fibroblasts", "Macrophages", "Neutrophils", "NK_cells", "Plasma_cells")
color_palette <- c("B_cells"="#e85c59", "CD4_CD8_T_cells"="#d99a00", "DC_cells"="#a9a9a9", "Endothelial_cells"="#a3b300", "Epithelial_cells"="#00b553", "Fibroblasts"="#00bfc4", "Macrophages"="#00b8e5", "Neutrophils"="#619cff", "NK_cells"="#d378f0", "Plasma_cells"="#f876aa")
obj_rna$celltype_plot <- factor(obj_rna$celltype_plot, levels = celltype_order)

# --- 步骤 4: 筛选并绘图 (与之前相同) ---
bhlhe40_expression <- GetAssayData(obj_rna, assay = "RNA", slot = "data")["BHLHE40", ]
cells_with_expression <- names(bhlhe40_expression[bhlhe40_expression > 0])
obj_rna_subset <- subset(obj_rna, cells = cells_with_expression)

vln_plot_filtered <- VlnPlot(
  object = obj_rna_subset,
  features = "BHLHE40",
  group.by = "celltype_plot", 
  pt.size = 0, 
  cols = color_palette 
) +
  labs(
    title = "BHLHE40 Expression (Cells with >0 Expression)",
    x = NULL, 
    y = "Expression Level"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))
print(vln_plot_filtered)
ggsave("BHLHE40_vlnplot.pdf", plot = vln_plot_filtered, width = 8, height = 6)

ggsave("BHLHE40_vlnplot.png", plot = vln_plot_filtered, width = 8, height = 6, dpi = 300)
