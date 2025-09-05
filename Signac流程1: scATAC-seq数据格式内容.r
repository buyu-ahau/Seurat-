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