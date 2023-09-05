# - Generate figures ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(lubridate)
source("scripts/figure_functions.R")

# make output directory
if (!dir.exists("figures/panels")) dir.create("figures/panels")

# - Read in general data ------------------------------------------------------
celltype_colors <- read.csv("data/celltype_colors.csv") %>% tibble::deframe()

# - Sun1 dataset --------------------------------------------------------------
# all cells
srt <- readRDS("results/sun1/all_cells.rds")

# figure S2B/C -- genes & UMIs per sample
p1 <- plot_metadata(srt, "nFeature_RNA", fills = "gray")
p2 <- plot_metadata(srt, "nCount_RNA", fills = "gray")

# figure S2D -- UMAP plot of celltype
p3 <- plot_umap(srt, groupby = "Celltype", colors = celltype_colors, legend = FALSE)

# figure S2E -- percentage of each celltype
p4 <- plot_percentage(srt, "Celltype", celltype_colors, flipped = TRUE)

# figure S2F -- Lepr expression across celltypes
p5 <- plot_gene_average(srt, "Lepr", groupby = "Celltype", slot = "data") +
  theme(axis.text.y = element_blank(),
        legend.key.width = unit(0.1, "in"))

# combine all all-cell plots
fig <- p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 5, widths = c(0.8, 0.8, 3, 1.2, 0.6))
ggsave("figures/panels/S2A-E.png", fig, width = 6.5, height = 3.2, units = "in", dpi = 400)

# neurons
srt <- readRDS("results/sun1/neurons.rds")
neuron_colors <- read.csv("results/sun1/neuron_colors.csv") %>% tibble::deframe()

p1 <- plot_umap(srt, label_clusters = FALSE, colors = neuron_colors)
p2 <- plot_gene_umap(srt, "Tcf7l2")
p3 <- plot_gene_umap(srt, "En1")
p4 <- plot_gene_umap(srt, "Adora2a")
fig <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4, widths = c(1, 1.1, 1.1, 1.1))
ggsave("figures/panels/S2F-I.png", fig, width = 7.5, height = 2, units = "in", dpi = 400)

# hypothalamic neurons
hypo <- readRDS("results/sun1/hypo.rds")
colors <- read.csv("results/sun1/hypo_cluster_colors.csv") %>% tibble::deframe()

# rename
new_names <- c(
  "Glp1r/Ebf1" = "Glp1r",
  "Agrp/Npy" = "Agrp",
  "Nts/Crh" = "Nts"
)
Idents(hypo) <- recode(Idents(hypo), !!!new_names)
names(colors) <- recode(names(colors), !!!new_names)

# figure 1C -- UMAP plot of Sun1 clusters
p1 <- plot_umap(hypo, color = colors)
ggsave("figures/panels/1C.png", p1, width = 3, height = 3, units = "in",
       dpi = 400, bg = "white")

# figure 1D -- marker genes
markers <- read.csv("results/sun1/markers_hypo.csv")
markers$cluster <- factor(markers$cluster, levels = levels(Idents(hypo)))
top10 <- markers %>% filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(gene) %>%
  slice(1) %>% ungroup() %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:10) %>%
  pull(gene)

p1 <- ScaleData(hypo, top10) %>%
  plot_heatmap(genes = top10, bars = TRUE, colors = colors)

# figures 1E and F
p2 <- plot_percentage(hypo, flipped = TRUE, fills = colors, reverse = TRUE,
                      transparent = c("GABA", "GLUT1", "GLUT2", "GLUT3")) +
  ylab(NULL)

p3 <- get_lepr(hypo, srt) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Cluster") %>%
  mutate(Cluster = factor(Cluster, levels = levels(Idents(hypo))[length(unique(Idents(hypo))):1])) %>%
  ggplot(aes(x = "", y = Cluster, fill = `.`)) +
  geom_tile(color = "black", size = edge_width) +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 1, size = 6),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.height = unit(0.02, "in"),
        legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
        scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                             breaks = seq(-2, 2), labels = seq(-2, 2))

fig <- p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(1, 0.5, 0.1))
ggsave("figures/panels/1D-F.png", fig, width = 5.3, height = 3, units = "in",
       dpi = 400, bg = "white")

# - TRAP enrichment -----------------------------------------------------------
markers <- read.csv("results/sun1/trap_enriched_markers.csv")
markers$cluster <- factor(markers$cluster, levels = levels(Idents(hypo)))

markers$cluster <- recode(markers$cluster, "Glp1r" = "Glp1r/Ebf1", "Agrp" = "Agrp/Npy",
                          "GABA1" = "GABA", "GABA2" = "Pvalb")

genes <- markers %>% filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(gene) %>%
  slice(1) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:10) %>%
  pull(gene)

hypo <- ScaleData(hypo, genes)
p1 <- plot_heatmap(hypo, genes = genes, flipped = TRUE, bars = TRUE, colors = colors)
ggsave("figures/panels/2B.png", p1, width = 3, height = 2.2, units = "in", dpi = 400)

# TRAP leptin treatment
de <- read.csv("results/GSE162603/de.csv")
de <- filter(de, gene_name %in% rownames(srt))

genes <- c("Stat3", "Asb4", "Atf3", "Timp1", "Serpina3i", "Prokr2", "Nlrc5",
           "Socs3", "Traf3ip3", "Serpina3n", "Irf9")

# 2C -- volcano of DE genes
p1 <- plot_volcano(de, label = genes) + xlab(NULL)

# 2D -- expression of highlighted genes
p2 <- plot_gene_average(hypo, genes, reverse = TRUE, normalize = TRUE, slot = "data") +
  theme(legend.key.width = unit(0.1, "in"),
        axis.text.y = element_text(size = 7))

fig <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 0.7))
ggsave("figures/panels/2C-D.png", fig, width = 3.4, height = 2.6, units = "in", dpi = 400)

fig <- FeaturePlot(hypo, "Leptin") +
  theme_void() +
  scale_color_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080") +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.height = unit(0.02, "in"),
        legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
ggsave("figures/panels/2E.png", fig, width = 2.2, height = 2.4, units = "in", dpi = 400)

# - Mouse, rat, and macaque snRNA-seq -----------------------------------------
# all cells
mouse <- readRDS("results/conservation/mouse.rds")
rat <- readRDS("results/conservation/rat.rds")
nhp <- readRDS("results/conservation/nhp.rds")

# plot metadata
p1 <- plot_metadata(mouse, "nFeature_RNA", fills = c("gray"))
p2 <- plot_metadata(rat, "nFeature_RNA", fills = c("#AF6458", "#526A83"))
p3 <- plot_metadata(nhp, "nFeature_RNA", fills = c("#AF6458", "#526A83"))
p4 <- plot_umap(mouse, groupby = "Sample", colors = "gray", label_clusters = FALSE) + ggtitle("Mouse")
p5 <- plot_umap(rat, groupby = "Sample", colors = c("#AF6458", "#526A83"), label_clusters = FALSE) + ggtitle("Rat")
p6 <- plot_umap(nhp, groupby = "Sample", colors = c("#AF6458", "#526A83"), label_clusters = FALSE) + ggtitle("Macaque")

fig <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 6, widths = c(1, 2, 2, 4, 4, 4))
ggsave("figures/panels/S_species_metadata.png", fig,
       width = 7.5, height = 2.5,
       units = "in", dpi = 400)

# umap plot of celltypes and Snap25 expression
fig <- cowplot::plot_grid(
  plot_umap(mouse, groupby = "Celltype", colors = celltype_colors) + ggtitle("Mouse"),
  plot_umap(rat, groupby = "Celltype", colors = celltype_colors) + ggtitle("Rat"),
  plot_umap(nhp, groupby = "Celltype", colors = celltype_colors) + ggtitle("Macaque"),
  plot_feature(mouse, "Snap25"),
  plot_feature(rat, "Snap25"),
  plot_feature(nhp, "SNAP25"),
  ncol = 3
)
ggsave("figures/panels/SUMAP_species_celltype.png", fig,
       width = 6.5, height = 4.5,
       units = "in", dpi = 400)

# table of celltypes
df <- bind_rows(
  as.data.frame(table(mouse$Celltype)/ncol(mouse)*100) %>% mutate(Species = "Mouse"),
  as.data.frame(table(rat$Celltype)/ncol(rat)*100) %>% mutate(Species = "Rat"),
  as.data.frame(table(nhp$Celltype)/ncol(nhp)*100) %>% mutate(Species = "Macaque")
) %>%
  mutate(Species = factor(Species, levels = c("Mouse", "Rat", "Macaque")))
celltype_levels <- df %>% group_by(Var1) %>%
  summarize("mean" = mean(Freq)) %>%
  arrange(mean) %>%
  filter(Var1 != "Doublets") %>%
  pull(Var1)
p1 <- df %>% mutate(Var1 = factor(Var1, levels = celltype_levels)) %>%
  filter(!is.na(Var1)) %>%
  tidyr::complete(Var1, Species, fill = list("Freq" = 0)) %>%
  ggplot(aes(x = Species, y = Var1, fill = Freq)) +
    geom_tile() +
    theme_void() +
    scale_fill_gradient(low = "#fbe6c5", high = "#70284a", name = "%") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = element_text(hjust = 1, size = 8),
          legend.key.width = unit(0.1, "in"),
          legend.key.height = unit(0.02, "in"),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 6),
          legend.title = element_blank())
fig <- p1 + p2 + plot_layout(ncol = 1)
ggsave("figures/panels/S_conserved_props+Lepr.png", fig,
       width = 1, height = 4.5,
       units = "in", dpi = 400)

# Lepr expression
df <- bind_rows(
  sapply(split(mouse[["RNA"]]@data["Lepr", ], mouse$Celltype), mean) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    mutate(Species = "Mouse"),
  sapply(split(rat[["RNA"]]@data["Lepr", ], rat$Celltype), mean) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    mutate(Species = "Rat"),
  sapply(split(nhp[["RNA"]]@data["LEPR", ], nhp$Celltype), mean) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    mutate(Species = "Macaque")
) %>%
  mutate(Species = factor(Species, levels = c("Mouse", "Rat", "Macaque")))
p2 <- df %>% mutate(Var1 = factor(Var1, levels = celltype_levels)) %>%
  filter(!is.na(Var1)) %>%
  tidyr::complete(Var1, Species, fill = list("." = 0)) %>%
  ggplot(aes(x = Species, y = Var1, fill = `.`)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient(low = "#fbe6c5", high = "#70284a",
                      name = expression(italic("Lepr"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(hjust = 1, size = 8),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.02, "in"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave("figures/panels/S_conserved_celltype_props.png", p,
       width = 1, height = 2.5,
       units = "in", dpi = 400)

# neurons alone
mouse <- readRDS("results/conservation/mouse_neurons.rds")
rat <- readRDS("results/conservation/rat_neurons.rds")
nhp <- readRDS("results/conservation/nhp_neurons.rds")

# integrated dataset
merged <- readRDS("results/conservation/integrated_neurons.rds")

# load colors
mouse_colors <- load_colors("mouse")
rat_colors <- load_colors("rat")
nhp_colors <- load_colors("nhp")

p1 <- plot_umap(mouse, colors = mouse_colors, label = FALSE)
p2 <- plot_umap(rat, reduction = "ortho.umap", colors = rat_colors, label = FALSE)
p3 <- plot_umap(nhp, reduction = "ortho.umap", colors = nhp_colors, label = FALSE)
fig <- p1 + p2 + p3
ggsave("figures/panels/3A.png", fig, width = 2.1, height = 0.7, units = "in", dpi = 400)

# plot species heatmap
p1 <- plot_species_heatmap("mouse")
ggsave("figures/panels/3B.png", p1, width = 1.5, height = 2.6, units = "in", dpi = 400, bg = "white")
p2 <- plot_species_heatmap("rat") + theme(legend.position = "bottom")
ggsave("figures/panels/3C.png", p2, width = 1.5, height = 2.6, units = "in", dpi = 400, bg = "white")
p3 <- plot_species_heatmap("nhp") + theme(legend.position = "bottom")
ggsave("figures/panels/3D.png", p3, width = 1.5, height = 2.6, units = "in", dpi = 400, bg = "white")

# UMAP combined
p1 <- plot_umap_combined("mouse")
p2 <- plot_umap_combined("rat")
p3 <- plot_umap_combined("nhp")
fig <- p1 + p2 + p3
ggsave("figures/panels/3F-H.png", fig, width = 5.7, height = 1.9, units = "in",
       dpi = 400, bg = "white")

# plot merged UMAP by species
species_colors <- c(
  colorRampPalette(c("white", "#672044"))(101)[75],
  colorRampPalette(c("white", "#123f5a"))(101)[75],
  colorRampPalette(c("white", "#3d5941"))(101)[75]
)
to_plot <- sapply(split(Cells(merged), merged$Species), sample, 4000) %>% sample()
p <- bind_cols(
  data.frame("Species" = merged$Species[to_plot],
  as.data.frame(merged[["umap"]]@cell.embeddings[to_plot, ]))
) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Species)) +
  geom_point(stroke = 0, size = 0.6) +
  scale_color_manual(values = species_colors, guide = "none") +
  theme_void()
ggsave("figures/panels/3I.png", p, width = 2.2, height = 2.2, units = "in",
       dpi = 400, bg = "white")

# plot UMAP of merged dataset split by species
clusters <- c("Glp1r", "Irx5", "Agrp", "Nts", "Pomc", "Prlh",
              "Foxb1", "Nr5a1", "Ghrh", "KNDy", "Pvalb")

merged_colors <- vector("character", length = length(levels(Idents(merged))))
names(merged_colors) <- levels(Idents(merged))
merged_colors[clusters] <- colors[clusters]
merged_colors[merged_colors == ""] <- paste0("gray", round(seq(50, 90, length.out = sum(merged_colors == "")), 0))
p1 <- plot_merged("Sun1")
p2 <- plot_merged("Mouse")
p3 <- plot_merged("Rat")
p4 <- plot_merged("Macaque")
fig <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave("figures/panels/S4A-D.png", fig, width = 5.5, height = 5.5, units = "in", dpi = 400)

p1 <- plot_correlation()
ggsave("figures/panels/S4E-G.png", p1, width = 2, height = 5.5, units = "in", dpi = 400)

# remove Sun1 sample from merged object
merged <- subset(merged, Sample != "Mouse_Sample_1776-AT-1")

genes <- c("Glp1r", "Irx5", "Npy", "Nts", "Pomc", "Prlh", "Foxb1", "Nr5a1", "Ghrh", "Kiss1", "Pvalb")
if (!all(genes %in% rownames(merged[["RNA"]]@scale.data))) {
  to_scale <- union(genes, rownames(merged[["RNA"]]@scale.data))
  merged <- ScaleData(merged, features = to_scale, split.by = "Sample")
}
p1 <- plot_species_genes(genes, clusters = clusters[length(clusters):1])

# plot Lepr expression
p2 <- read.csv("results/conservation/lepr_expression.csv") %>%
  mutate("X" = factor(X, levels = X[nrow(.):1])) %>%
  tidyr::pivot_longer(-X, names_to = "Species", values_to = "Lepr") %>%
  mutate(Species = factor(Species, levels = c("Mouse", "Rat", "Macaque"))) %>%
  ggplot(aes(x = Species, y = X, fill = Lepr)) +
  geom_tile() +
  scale_fill_gradient(low = "#fbe6c5", high = "#70284a",
                      name = expression(italic("Lepr"))) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(hjust = 1, size = 8),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.02, "in"),
        legend.key.width = unit(0.1, "in"))

fig <- p1 + p2 + plot_layout(ncol = 2, widths = c(4, 0.8))
ggsave("figures/panels/3J-K.png", fig, width = 5, height = 2.5, units = "in", dpi = 400)


# - Glp1rCre;LeprFlox ---------------------------------------------------------
plot_coexpression <- function(object, gene1, gene2) {
  xlab <- paste0(gene1, "/", gene2)
  get_coexpression(object, gene1, gene2) %>%
    mutate("Var2" = factor(Var2, levels = levels(Idents(object))[length(unique(Idents(object))):1])) %>%
    ggplot(aes(x = Freq, y = Var2, fill = Var2)) +
    geom_col(color = "black", size = edge_width) +
    theme_classic() +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    xlab(paste0("<i>", xlab, "</i> cells")) +
    ylab(NULL) +
    theme(axis.text.x = element_text(size = 7, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.x = ggtext::element_markdown(size = 8))
}
#ggsave("figures/panels/4A.png", p1, width = 1.8, height = 2.2, units = "in", dpi = 400)

# new figures
p <- cowplot::plot_grid(
  plot_coexpression(hypo, "Glp1r", "Lepr"),
  plot_coexpression(hypo, "Agrp", "Lepr"),
  plot_coexpression(hypo, "Pomc", "Lepr"),
  ncol = 3
)
ggsave("figures/panels/counts.png", p, width = 1.8*3, height = 2.2, units = "in", dpi = 400)

# alternative
p <- cowplot::plot_grid(
  plot_coexpression(hypo, "Glp1r", "Lepr"),
  plot_coexpression(hypo, "Npy", "Lepr"),
  plot_coexpression(hypo, "Pomc", "Lepr"),
  ncol = 3
)
ggsave("figures/panels/counts_alt.png", p, width = 1.8*3, height = 2.2, units = "in", dpi = 400)