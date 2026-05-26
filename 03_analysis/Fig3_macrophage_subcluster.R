library(future)
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(ggprism)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scales)

plan("multicore", workers = 10)
options(future.globals.maxSize = 1024^3 * 200)

set.seed(123)

theme_set(ggpubr::theme_pubr() + theme(legend.position = "top"))

dir.create("output", showWarnings = FALSE)

colors1 <- c(
  "#2B3D26", "#F3C300", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
  "#008856", "#E68FAC", "#AAF400", "#00A5FF", "#604E97", "#FFE09D", "#B3446C", "#1CFFCE",
  "#882D17", "#8DB600", "#654522", "#E25822", "#F6222E", "#FE00FA",
  "#F6A600", "#3283FE", "#FEAF16", "#B00068", "#DCD300", "#90AD1C", "#2ED9FF", "#DEA0FD",
  "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D", "#B10DA1", "#FBE426",
  "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", "#0067A5", "#BDCDFF",
  "#822E1C", "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB",
  "#D100E5", "#60E1E0", "#F66196", "#7870B3", "#84DE02", "#36BDA3", "#875692", "#3D52D5",
  "#B4A8BD", "#FDDC22", "#16FF32", "#1C9963", "#6ECD6E", "#6A7F7A", "#0A437A", "#FFA0F2",
  "#CCAA35", "#F99379", "#D300E7", "#9003C7", "#FF5733", "#0075DC", "#9EFD38", "#D8A903",
  "#00A68C", "#EF798A", "#F0E68C", "#15D3A1", "#FF6E4A", "#1E3A4E", "#0D98BA", "#2E5894",
  "#9C2542", "#FF8D00", "#69359C", "#C62D42", "#ADD8E6", "#E6E200", "#66FF00", "#BF94E4",
  "#00FF6F", "#EE34D2", "#BBFFFF", "#FFD700", "#C2B280", "#8DB600", "#654522", "#E25822",
  "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068"
)

dim_usage  <- 20
resolution_macro  <- 0.65
resolution_inflam <- 0.5
pt_size    <- 0.3
logFC_cut  <- 0.58
pval_cut   <- 0.05

feature_marker <- list(gy = c("newRNA_ratio"))

make_colors <- function(obj, group) {
  colorRampPalette(brewer.pal(12, "Paired"))(length(unique(obj@meta.data[[group]])))
}

plot_features <- function(obj, markers, ncol = 1, min.cutoff = NA, max.cutoff = NA) {
  p1 <- FeaturePlot(obj, features = unlist(markers), order = TRUE,
                    combine = FALSE, reduction = "umap",
                    min.cutoff = min.cutoff, max.cutoff = max.cutoff)
  fix.sc <- scale_color_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
  p2 <- lapply(p1, function(x) x + fix.sc)
  plot_grid(plotlist = p2, ncol = ncol)
}

all_obj    <- readRDS("all.rds")
Macrophage <- subset(all_obj, celltype1 == "Macrophage")
DefaultAssay(Macrophage) <- "RNA"
table(Macrophage$seurat_clusters)

Macrophage <- RunHarmony(Macrophage, "orig.ident")
Macrophage <- RunUMAP(Macrophage, reduction = "harmony", dims = 1:dim_usage)
Macrophage <- FindNeighbors(Macrophage, reduction = "harmony", dims = 1:dim_usage)
Macrophage <- FindClusters(Macrophage, resolution = resolution_macro)

options(repr.plot.height = 5, repr.plot.width = 5)
DimPlot(Macrophage, reduction = "umap", label = TRUE,
        cols = make_colors(Macrophage, "seurat_clusters"),
        pt.size = pt_size, group.by = "seurat_clusters")

new.cluster.ids <- c("Inflam_TAMs", "Reg_TAMs", "Prolif_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs",
                     "Reg_TAMs", "DC_cells", "Reg_TAMs", "INF_TAMs", "DC_cells")
Idents(Macrophage) <- Macrophage$seurat_clusters
names(new.cluster.ids) <- levels(Macrophage)
Macrophage <- RenameIdents(Macrophage, new.cluster.ids)
Macrophage$celltype1 <- Idents(Macrophage)

custom_color <- make_colors(Macrophage, "celltype1")
DimPlot(Macrophage, reduction = "umap", label = TRUE,
        cols = custom_color, pt.size = pt_size, group.by = "celltype1")

print(plot_features(Macrophage, feature_marker))

Inflam_TAMs <- subset(Macrophage, celltype1 == "Inflam_TAMs")
DefaultAssay(Inflam_TAMs) <- "RNA"
table(Inflam_TAMs$seurat_clusters)

Inflam_TAMs <- RunHarmony(Inflam_TAMs, "orig.ident")
Inflam_TAMs <- RunUMAP(Inflam_TAMs, reduction = "harmony", dims = 1:dim_usage)
Inflam_TAMs <- FindNeighbors(Inflam_TAMs, reduction = "harmony", dims = 1:dim_usage)
Inflam_TAMs <- FindClusters(Inflam_TAMs, resolution = resolution_inflam)

options(repr.plot.height = 5, repr.plot.width = 5)
DimPlot(Inflam_TAMs, reduction = "umap", label = TRUE,
        cols = make_colors(Inflam_TAMs, "seurat_clusters"),
        pt.size = pt_size, group.by = "seurat_clusters")

print(plot_features(Inflam_TAMs, feature_marker, min.cutoff = 0, max.cutoff = 0.5))

new.cluster.ids <- c("L-Inflam_TAMs", "H-Inflam_TAMs", "H-Inflam_TAMs",
                     "H-Inflam_TAMs", "H-Inflam_TAMs", "H-Inflam_TAMs")
Idents(Inflam_TAMs) <- Inflam_TAMs$seurat_clusters
names(new.cluster.ids) <- levels(Inflam_TAMs)
Inflam_TAMs <- RenameIdents(Inflam_TAMs, new.cluster.ids)
Inflam_TAMs$celltype1 <- Idents(Inflam_TAMs)

custom_color <- make_colors(Inflam_TAMs, "celltype1")
DimPlot(Inflam_TAMs, reduction = "umap", label = TRUE,
        cols = custom_color, pt.size = pt_size, group.by = "celltype1")

df <- Inflam_TAMs@meta.data
df$rate_percent <- df$newRNA_ratio * 100
ggplot(df, aes(x = celltype1, y = rate_percent, fill = celltype1)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_color) +
  theme_classic() +
  labs(x = "Cell Type", y = "New-to-total RNA ratios") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 50))

table(Inflam_TAMs$time)
prop.table(table(Idents(Inflam_TAMs)))
table(Idents(Inflam_TAMs), Inflam_TAMs$time)

Cellratio <- prop.table(table(Inflam_TAMs$celltype1, Inflam_TAMs$time), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c("NC", "D3", "D5", "D7", "D9", "D11", "D20"))

ggplot(Cellratio, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_stratum(width = 0.7, color = "white") +
  geom_alluvium(alpha = 0.5, width = 0.7, color = "white", size = 1, curve_type = "linear") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Samples", y = "Relative Abundance(%)", fill = "group") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  theme_prism(palette = "candy_bright", base_fontface = "plain", base_family = "serif",
              base_size = 16, base_line_size = 0.8, axis_text_angle = 45) +
  scale_fill_manual(values = custom_color) +
  theme(legend.position = "top")

highvslow <- FindMarkers(object = Inflam_TAMs, ident.1 = "H-Inflam_TAMs", ident.2 = "L-Inflam_TAMs",
                         min.pct = 0.25)
highvslow$symbol <- rownames(highvslow)

k1 <- (highvslow$p_val_adj < pval_cut) & (highvslow$avg_log2FC < -logFC_cut)
k2 <- (highvslow$p_val_adj < pval_cut) & (highvslow$avg_log2FC >  logFC_cut)
highvslow <- mutate(highvslow, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
highvslow <- subset(highvslow, p_val_adj != 0)

p_volcano <- ggplot(data = highvslow, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.4, size = 3.5, aes(color = change)) +
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("blue4", "grey", "red3")) +
  geom_vline(xintercept = c(-logFC_cut, logFC_cut), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(pval_cut),          lty = 4, col = "black", lwd = 0.8) +
  theme_classic()

up_data   <- filter(highvslow, change == "up")   %>% distinct(symbol, .keep_all = TRUE) %>% top_n(5, -log10(p_val_adj))
down_data <- filter(highvslow, change == "down")  %>% distinct(symbol, .keep_all = TRUE) %>% top_n(5, -log10(p_val_adj))

p_volcano_labeled <- p_volcano +
  geom_text_repel(data = up_data,   aes(x = avg_log2FC, y = -log10(p_val_adj), label = symbol)) +
  geom_text_repel(data = down_data, aes(x = avg_log2FC, y = -log10(p_val_adj), label = symbol))

p_volcano_labeled

gene_df <- bitr(subset(highvslow, change == "up")$symbol,
                fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ego_BP <- enrichGO(gene_df$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
                   pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05, readable = TRUE)

barplot(ego_BP, showCategory = 10)
