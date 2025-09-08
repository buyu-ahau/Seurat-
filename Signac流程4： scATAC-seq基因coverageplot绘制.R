# ---适配GTF注释以匹配Seurat对象并绘图 ---
# 目的:
# 1. 遵循最正确的原则：保持原始Seurat对象的坐标系统不变。
# 2. 将从外部GTF文件加载的基因注释，强制转换为与Seurat对象一致的Ensembl命名风格 ('1', 23456'13'...)。
# 3. 最终成功为目标基因BHLHE40绘制覆盖度图。
# ----------------------------------------------------------------
library(Signac)
library(Seurat)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
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
cat("--- 步骤 4: 正在为基因BHLHE40绘制覆盖度图...\n")
# 现在，所有坐标系统都已统一，绘图将成功执行
p <- CoveragePlot(
  object = obj_relinked,
  region = "BHLHE40",
  group.by = "seurat_clusters",
  extend.upstream = 2000,
  extend.downstream = 1000
)
# 美化并展示图像
p_styled <- p + theme_classic(base_size = 12) +
            labs(title = "BHLHE40 Gene Locus Chromatin Accessibility",
                 subtitle = "Grouped by Cell Type (seurat_clusters)")
print(p_styled)
# Save as PDF
ggsave(filename = "BHLHE40_coverage_plot.pdf", plot = p, width = 12, height = 10)
# Save as PNG with high resolution
ggsave(filename = "BHLHE40_coverage_plot.png", plot = p, width = 12, height = 10, dpi = 300)