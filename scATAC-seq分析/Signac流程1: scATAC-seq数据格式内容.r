
> file_path <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/daa52225-obj.Rda"
> load(file_path)
> object_names <- ls()
print("已加载的对象名称:")
print(object_names)
[1] "已加载的对象名称:"
[1] "file_path" "obj"      
--- 1. 对象的类别 ---
[1] "Seurat"
attr(,"package")
[1] "SeuratObject"

--- 2. 对象摘要信息 ---
Loading required namespace: SeuratObject
Loading required namespace: Signac
An object of class Seurat 
281627 features across 34310 samples within 2 assays 
Active assay: ATAC (259795 features, 64998 variable features)
 2 layers present: counts, data
 1 other assay present: ACTIVITY
 9 dimensional reductions calculated: lsi, tsne, tsne_ATAC, umap, umap_ATAC, harmony_ATAC, harmony, tsne_harmony_ATAC, umap_harmony_ATAC

--- 3. 对象的内部结构 ---
Formal class 'Seurat' [package "SeuratObject"] with 13 slots
  ..@ assays      :List of 2
  ..@ meta.data   :'data.frame':        34310 obs. of  11 variables:
  ..@ active.assay: chr "ATAC"
  ..@ active.ident: Factor w/ 9 levels "B_cells","CD4_CD8_T_cells",..: 6 2 9 2 6 9 9 8 2 6 ...
  .. ..- attr(*, "names")= chr [1:34310] "Mhp_1_AAACGAAAGGTCACTT-1" "Mhp_1_AAACGAAAGTCGACCC-1" "Mhp_1_AAACGAAAGTCGTACT-1" "Mhp_1_AAACGAACAATTCAGC-1" ...
  ..@ graphs      :List of 2
  ..@ neighbors   : list()
  ..@ reductions  :List of 9
  ..@ images      : list()
  ..@ project.name: chr "SeuratProject"
  ..@ misc        :List of 5
  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
  ..@ commands    :List of 6
  ..@ tools       : list()

--- 4. Seurat/Signac 对象特有信息 ---

   a. 包含的 Assays:
[1] "ATAC"     "ACTIVITY"

   b. 元数据预览 (前6个细胞):
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

   c. 降维结果:
[1] "lsi"               "tsne"              "tsne_ATAC"        
[4] "umap"              "umap_ATAC"         "harmony_ATAC"     
[7] "harmony"           "tsne_harmony_ATAC" "umap_harmony_ATAC"