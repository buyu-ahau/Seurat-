load("daa52225-obj.Rda")
library(Seurat)
library(ggplot2)
library(dplyr)

plotting_df <- obj@reductions$umap_harmony_ATAC@cell.embeddings %>%
  as.data.frame() %>%
  select(UMAP_1 = 1, UMAP_2 = 2) %>%
  merge(obj@meta.data, by = 0) %>%
  rename(Cluster = seurat_clusters, Barcode = Row.names) %>%
  select(UMAP_1, UMAP_2, Cluster)

color_palette <- c(
  "#919ac2", "#ffac98", "#70a4c8", "#a5a9af", "#63917d",
  "#dbd1b4", "#6e729a", "#9ba4bd", "#c5ae5f", "#b9b8d6"
)

umap_plot <- ggplot(plotting_df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.2, shape = 16) +
  scale_color_manual(values = color_palette) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5,
      arrow = arrow(length = unit(0.2, "cm"))
    ),
    aspect.ratio = 1,
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# 保存图片
ggsave(
  filename = "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/scATAC-seq云雾图UMAP.png",
  plot = umap_plot,
  height = 6,
  width = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "scATAC-seq云雾图UMAP.pdf",
  plot = umap_plot,
  height = 6,
  width = 7
)





--- 1. 宏观概览 ---
Loading required package: Signac
An object of class Seurat 
281627 features across 34310 samples within 2 assays 
Active assay: ATAC (259795 features, 64998 variable features)
 2 layers present: counts, data
 1 other assay present: ACTIVITY
 9 dimensional reductions calculated: lsi, tsne, tsne_ATAC, umap, umap_ATAC, harmony_ATAC, harmony, tsne_harmony_ATAC, umap_harmony_ATAC


--- 2. 数据层 (Assays) 详细信息 ---
a. 对象中包含的Assay(s):
[1] "ATAC"     "ACTIVITY"

b. 当前活动的Assay:
[1] "ATAC"

c. 活动Assay ('ATAC') 详情:
ChromatinAssay data with 259795 features for 34310 cells
Variable features: 64998 
Genome: 
Annotation present: TRUE 
Motifs present: FALSE 
Fragment files: 6 


--- 3. 细胞元数据 (Metadata) 分析 ---
a. 元数据包含的列名:
 [1] "orig.ident"            "nCount_ATAC"           "nFeature_ATAC"        
 [4] "is__cell_barcode"      "peak_region_fragments" "passed_filters"       
 [7] "nCount_ACTIVITY"       "nFeature_ACTIVITY"     "ATAC_snn_res.0.5"     
[10] "seurat_clusters"       "beforeInteg.cluster"  

b. 元数据内容预览 (前6个细胞):
                         orig.ident nCount_ATAC nFeature_ATAC is__cell_barcode
Mhp_1_AAACGAAAGGTCACTT-1      Mhp_1        3110          2825                1
Mhp_1_AAACGAAAGTCGACCC-1      Mhp_1        6139          4896                1
Mhp_1_AAACGAAAGTCGTACT-1      Mhp_1        2496          2424                1
Mhp_1_AAACGAACAATTCAGC-1      Mhp_1        2083          1692                1
Mhp_1_AAACGAACACTAGGAG-1      Mhp_1        1495          1459                1
Mhp_1_AAACGAACAGAACTTC-1      Mhp_1        2167          1823                1
                         peak_region_fragments passed_filters nCount_ACTIVITY
Mhp_1_AAACGAAAGGTCACTT-1                  2246          10561            4936
Mhp_1_AAACGAAAGTCGACCC-1                  4992          11301            6575
Mhp_1_AAACGAAAGTCGTACT-1                  1907           6109            3265
Mhp_1_AAACGAACAATTCAGC-1                  1396           8275            3698
Mhp_1_AAACGAACACTAGGAG-1                  1039           4828            2271
Mhp_1_AAACGAACAGAACTTC-1                  1566           6873            3407
                         nFeature_ACTIVITY ATAC_snn_res.0.5 seurat_clusters
Mhp_1_AAACGAAAGGTCACTT-1              3214                0     Macrophages
Mhp_1_AAACGAAAGTCGACCC-1              3286                6 CD4_CD8_T_cells
Mhp_1_AAACGAAAGTCGTACT-1              2655                4    Plasma_cells
Mhp_1_AAACGAACAATTCAGC-1              2046                7 CD4_CD8_T_cells
Mhp_1_AAACGAACACTAGGAG-1              1981                0     Macrophages
Mhp_1_AAACGAACAGAACTTC-1              2051                4    Plasma_cells
                         beforeInteg.cluster
Mhp_1_AAACGAAAGGTCACTT-1                   0
Mhp_1_AAACGAAAGTCGACCC-1                   8
Mhp_1_AAACGAAAGTCGTACT-1                   0
Mhp_1_AAACGAACAATTCAGC-1                   7
Mhp_1_AAACGAACACTAGGAG-1                   0
Mhp_1_AAACGAACAGAACTTC-1                   6

c. 各细胞群细胞数量统计:

          B_cells   CD4_CD8_T_cells Endothelial_cells  Epithelial_cells 
             2901              9495              2473              1937 
      Fibroblasts       Macrophages       Neutrophils          NK_cells 
             1891              9344              1296              2697 
     Plasma_cells 
             2276 


--- 4. 降维结果 (Reductions) 列表 ---
a. 已计算的降维方法:
[1] "lsi"               "tsne"              "tsne_ATAC"        
[4] "umap"              "umap_ATAC"         "harmony_ATAC"     
[7] "harmony"           "tsne_harmony_ATAC" "umap_harmony_ATAC"

b. 降维结果可视化预览 (UMAP):
UMAP图已生成 (请在Plots窗口查看)。

--- 5. 基因组注释信息检查 ---
对象已包含基因组注释信息。
GRanges object with 6 ranges and 21 metadata columns:
      seqnames              ranges strand |   source        type     score
         <Rle>           <IRanges>  <Rle> | <factor>    <factor> <numeric>
  [1]        1 226161299-226217308      - |  ensembl gene               NA
  [2]        1 226188008-226217308      - |  ensembl transcript         NA
  [3]        1 226217188-226217308      - |  ensembl exon               NA
  [4]        1 226217188-226217253      - |  ensembl CDS                NA
  [5]        1 226217251-226217253      - |  ensembl start_codon        NA
  [6]        1 226202702-226202806      - |  ensembl exon               NA
          phase            gene_id gene_version   gene_name gene_source
      <integer>        <character>  <character> <character> <character>
  [1]      <NA> ENSSSCG00000028996            4     ALDH1A1     ensembl
  [2]      <NA> ENSSSCG00000028996            4     ALDH1A1     ensembl
  [3]      <NA> ENSSSCG00000028996            4     ALDH1A1     ensembl
  [4]         0 ENSSSCG00000028996            4     ALDH1A1     ensembl
  [5]         0 ENSSSCG00000028996            4     ALDH1A1     ensembl
  [6]      <NA> ENSSSCG00000028996            4     ALDH1A1     ensembl
        gene_biotype      transcript_id transcript_version transcript_name
         <character>        <character>        <character>     <character>
  [1] protein_coding               <NA>               <NA>            <NA>
  [2] protein_coding ENSSSCT00000103363                  1     ALDH1A1-201
  [3] protein_coding ENSSSCT00000103363                  1     ALDH1A1-201
  [4] protein_coding ENSSSCT00000103363                  1     ALDH1A1-201
  [5] protein_coding ENSSSCT00000103363                  1     ALDH1A1-201
  [6] protein_coding ENSSSCT00000103363                  1     ALDH1A1-201
      transcript_source transcript_biotype exon_number            exon_id
            <character>        <character> <character>        <character>
  [1]              <NA>               <NA>        <NA>               <NA>
  [2]           ensembl     protein_coding        <NA>               <NA>
  [3]           ensembl     protein_coding           1 ENSSSCE00000336723
  [4]           ensembl     protein_coding           1               <NA>
  [5]           ensembl     protein_coding           1               <NA>
  [6]           ensembl     protein_coding           2 ENSSSCE00000191205
      exon_version         protein_id protein_version         tag
       <character>        <character>     <character> <character>
  [1]         <NA>               <NA>            <NA>        <NA>
  [2]         <NA>               <NA>            <NA>        <NA>
  [3]            2               <NA>            <NA>        <NA>
  [4]         <NA> ENSSSCP00000079993               1        <NA>
  [5]         <NA>               <NA>            <NA>        <NA>
  [6]            2               <NA>            <NA>        <NA>
      projection_parent_transcript
                       <character>
  [1]                         <NA>
  [2]                         <NA>
  [3]                         <NA>
  [4]                         <NA>
  [5]                         <NA>
  [6]                         <NA>
  -------
  seqinfo: 303 sequences from an unspecified genome; no seqlengths