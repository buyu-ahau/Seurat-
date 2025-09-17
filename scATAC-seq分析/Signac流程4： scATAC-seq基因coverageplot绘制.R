# ---适配GTF注释以匹配Seurat对象并绘图 ---
# 目的:
# 1. 遵循最正确的原则：保持原始Seurat对象的坐标系统不变。
# 2. 将从外部GTF文件加载的基因注释，强制转换为与Seurat对象一致的Ensembl命名风格 ('1', 23456'13'...)。
# 3. 最终成功为目标基因BHLHE40绘制覆盖度图。
# ----------------------------------------------------------------
###!!!!!!!!!conda activate signac_motif
library(Signac)
library(Seurat)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(qs)
# 加载Seurat对象
load("/disk192/users_dir/buyu/1/布宇/3/布宇scATAC-seq/2025.8.14/daa52225-obj-relinked-final.Rda")
# 从本地GTF文件导入基因注释
gtf_path <- "/disk192/users_dir/buyu/2/参考基因组/Sus_scrofa.Sscrofa11.1.114.gtf.gz"
pig_annotation_from_gtf <- rtracklayer::import(gtf_path)
# --- 步骤 2: 【关键修正】将GTF注释的风格适配为Seurat对象的风格 ---
# 首先，我们确认一下Seurat对象内部的命名风格 
print(head(seqlevels(granges(obj_relinked[['ATAC']]))))
# 然后，将GTF注释的命名风格强制转换为Ensembl风格，以匹配Seurat对象
# 这个命令会移除所有'chr'前缀，确保两者一致
seqlevelsStyle(pig_annotation_from_gtf) <- "Ensembl"
print(head(seqlevels(pig_annotation_from_gtf)))
# --- 步骤 3: 将适配好的注释更新到Seurat对象 ---
# 此时，obj_relinked对象本身没有被修改，我们只是给它加上了格式完全匹配的注释
Annotation(obj_relinked) <- pig_annotation_from_gtf
# --- 步骤 4: 绘制BHLHE40的Coverage Plot ---
# ===================================================================
# 步骤 1: 定义文件路径并读取您的数据
# ===================================================================
# 这是您提供的文件路径，请再次确认它是否准确无误
file_path <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/obj_relinked_annotated.rds"
# 我们将加载的对象命名为 obj_atac
print("正在加载数据，请稍候...")
obj_atac <- readRDS(file_path)
# (推荐) 打印对象的基本信息，确认数据已成功加载
print("数据加载成功！对象摘要如下：")
print(obj_atac)

print(Annotation(obj_atac))
# 将基因名称 "BHLHE40" 赋值给一个变量 roi
roi <- "BHLHE40"
print(paste("已设定目标基因为:", roi))
# 调用 CoveragePlot 函数进行绘图
print("正在生成 Coverage Plot...")
p_extended <- CoveragePlot(
  object = obj_atac,
  region = roi,
  extend.upstream = 2000,   # <--- 新增：向上游（左侧）扩展 2000 bp
  extend.downstream = 2000, # <--- 新增：向下游（右侧）扩展 2000 bp
  annotation = TRUE,        # 显示基因结构注释
  peaks = TRUE,             # 显示已识别的 Peaks
  tile = TRUE               # 显示 Tile 图，展示片段密度
)
print(p_extended)
ggsave("BHLHE40_CoveragePlot.png", plot = p_extended, width = 8, height = 6, dpi = 300)
ggsave("BHLHE40_CoveragePlot.pdf", plot = p_extended, width = 8, height = 6)





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