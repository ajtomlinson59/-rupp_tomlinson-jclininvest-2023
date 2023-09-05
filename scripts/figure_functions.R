# - Functions to plot figures -------------------------------------------------

library(ggplot2)
library(stringr)
library(dplyr)
library(Seurat)
library(ggrepel)
options(bitmapType = "cairo")

# col & tile edge width
edge_width <- 0.5

#' Plot given metadata as a violin
plot_metadata <- function(object, metadata, groupby = "Sample", fills = NULL) {
  if (metadata == "nFeature_RNA") title = "Genes"
  if (metadata == "nCount_RNA") title = "UMIs"
  p <- VlnPlot(object, metadata, group.by = groupby, pt.size = 0) +
    scale_y_continuous(trans = "log2", limits = c(1, NA)) +
    xlab(NULL) +
    ggtitle(title) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6, color = "black"),
          plot.title = element_text(size = 8, hjust = 0.5))
  if (!is.null(fills)) p <- p + scale_fill_manual(values = fills)
  return(p)
}

#' Plot cells in UMAP space based on label in `groupby`
plot_umap <- function(object, reduction = "umap", colors = NULL, groupby = NULL,
                      legend = FALSE, label_clusters = TRUE) {
  df <- bind_cols(
    as.data.frame(object[[reduction]]@cell.embeddings[, 1:2]),
    data.frame("Cluster" = Idents(object))
  )
  colnames(df)[1:2] <- c("Dim1", "Dim2")
  # cluster centers
  if (label_clusters) {
    centers <- df %>% group_by(Cluster) %>%
      summarize("Dim1" = median(Dim1), "Dim2" = median(Dim2))
  }
  # plot
  p <- ggplot(df, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = Cluster), stroke = 0, size = 0.6) +
    theme_void() +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)
  if (!legend) p <- p + theme(legend.position = "none")
  if (label_clusters) {
    p <- p + geom_text(data = centers, aes(label = Cluster), size = 2.5)
  }
  p
}

#' Plot percentage of a given object
plot_percentage <- function(object, groupby = "Clusters", fills = NULL,
                            flipped = FALSE, reverse = FALSE, transparent = NULL) {
  if (groupby == "Clusters") {
    df <- table(Idents(object), object$Sample)
  } else {
    df <- table(object@meta.data[, groupby], object$Sample)
  }
  # calculate percentage for each sample and format data for plotting
  df <- df %>%
    as.data.frame() %>%
    group_by(Var2) %>%
    mutate("Percent" = Freq / sum(Freq) * 100)
  if (length(unique(df$Var2)) > 1) {
    df <- df %>%
      group_by(Var1) %>%
      summarize("Avg" = mean(Percent), "SEM" = sd(Percent)/sqrt(n()))
  } else {
    df <- rename(df, "Avg" = Percent) %>% mutate("SEM" = NA)
  }
  if (reverse) df$Var1 <- factor(df$Var1, levels = levels(df$Var1)[length(levels(df$Var1)):1])
  # plot
  p <- ggplot(df, aes(x = Var1, y = Avg)) +
    geom_col(aes(fill = Var1, alpha = Var1 %in% transparent),
             color = "black", size = edge_width) +
    geom_errorbar(aes(ymin = Avg - SEM, ymax = Avg + SEM), width = 0.4) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(NULL) +
    ylab("%") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 6),
          axis.title = element_text(size = 8))
  if (flipped) p <- p + coord_flip()
  if (!is.null(fills)) p <- p + scale_fill_manual(values = fills)
  if (!is.null(transparent)) {
    p <- p + scale_alpha_manual(values = c(1, 0.5))
  } else {
    p <- p + scale_alpha_manual(values = 1)
  }
  p
}

#' Get gene average across clusters (or other grouping)
get_gene_average <- function(object, gene, groupby = NULL, slot = "scale.data") {
  if (!gene %in% rownames(slot(object[["RNA"]], slot))) {
    stop(gene, " is not in the ", slot, " slot")
  }
  expr <- slot(object[["RNA"]], slot) %>% .[gene, ]
  if (is.null(groupby)) {
    data <- split(expr, Idents(object))
  } else {
    data <- split(expr, object@meta.data[, groupby])
  }
  sapply(data, mean)
}

#' Plot gene averages
plot_gene_average <- function(object, genes, normalize = FALSE, reverse = FALSE, ...) {
  mtx <- sapply(genes, function(x) get_gene_average(object, x, ...))
  if (normalize) mtx <- apply(mtx, 2, function(x) x / max(x))
  df <- mtx %>% as.data.frame() %>%
    tibble::rownames_to_column("Cluster") %>%
    tidyr::pivot_longer(-Cluster, names_to = "Gene", values_to = "Expr") %>%
    mutate(Gene = factor(Gene, levels = genes))
  if (reverse) {
    df <- mutate(df, Cluster = factor(Cluster, levels = rownames(mtx)[nrow(mtx):1]))
  } else {
    df <- mutate(df, Cluster = factor(Cluster, levels = rownames(mtx)))
  }
  p <- ggplot(df, aes(x = Gene, y = Cluster, fill = Expr)) +
    geom_tile(color = "black", size = edge_width) +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     face = "italic", size = 7),
          axis.text.y = element_text(hjust = 1, size = 8),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.height = unit(0.02, "in"),
          legend.key.width = unit(0.5, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank())
  if (any(mtx < 0)) {
    p <- p + scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                                  breaks = seq(-2, 2), labels = seq(-2, 2))
  } else {
    p <- p + scale_fill_gradient(low = "#f6edbd", high = "#008080",
                                 breaks = c(0, 1), labels = c(0, 1))
  }
  p
}

#' Plot a feature in UMAP space
plot_feature <- function(object, feature, reduction = "umap", slot = "data") {
  if (!feature %in% rownames(slot(object[["RNA"]], slot))) {
    stop(feature, " not in ", slot, " slot")
  }
  df <- bind_cols(
    as.data.frame(object[[reduction]]@cell.embeddings[, 1:2]),
    data.frame("Expr" = slot(object[["RNA"]], slot)[feature, ])
  )
  colnames(df)[1:2] <- c("Dim1", "Dim2")
  # plot
  p <- ggplot(df, aes(x = Dim1, y = Dim2, color = Expr)) +
    geom_point(stroke = 0, size = 0.6) +
    theme_void() +
    ggtitle(feature) +
    scale_color_gradient(low = "#fbe6c5", high = "#70284a") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 8),
          legend.key.width = unit(0.02, "in"),
          legend.key.height = unit(0.1, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank())
  p
}

#' Plot gene expression across cells
plot_heatmap <- function(object, genes, limit = 2, slot = "scale.data",
                         n_cells = 1000, flipped = FALSE, bars = TRUE,
                         colors = NULL) {
  mtx <- slot(object[["RNA"]], slot)
  if (!all(genes %in% rownames(mtx))) {
    stop("Not all genes are in the ", slot, " slot")
  }
  if (ncol(mtx) > n_cells) mtx <- mtx[, sample(colnames(mtx), n_cells)]
  cell_order <- Idents(object)[colnames(mtx)] %>% sort() %>% names()
  if (any(abs(mtx) > limit)) {
    mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -limit, -limit, ifelse(x > limit, limit, x)))
  }
  # make tidy
  df <- as.data.frame(mtx) %>% tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expr")
  if (flipped) {
    df <- df %>%
      mutate(Gene = factor(Gene, levels = genes[length(genes):1])) %>%
      mutate(Cell = factor(Cell, levels = cell_order))
  } else {
    df <- df %>%
      mutate(Gene = factor(Gene, levels = genes)) %>%
      mutate(Cell = factor(Cell, levels = cell_order[length(cell_order):1]))
  }
  # plot
  p <- ggplot(df, aes(x = Gene, y = Cell, fill = Expr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                                  breaks = seq(-2, 2), labels = seq(-2, 2)) +
    theme_void() +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.key.height = unit(0.02, "in"),
          legend.key.width = unit(0.1, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank())
  if (bars) {
    if (is.null(colors)) {
      colors <- sample(grDevices::colors(), length(unique(Idents(object))))
      names(colors) <- levels(Idents(object))
    }
    height <- length(unique(df$Gene)) * 0.03
    for (i in unique(Idents(object))) {
      p <- p + annotate(
        "rect",
        xmin = 0.5 - height,
        xmax = 0.5,
        ymin = min(which(levels(df$Cell) %in% Cells(object)[Idents(object) == i])),
        ymax = max(which(levels(df$Cell) %in% Cells(object)[Idents(object) == i])),
        fill = colors[i],
        color = "black",
        size = edge_width
      )
    }
  }
  if (flipped) p <- p + coord_flip()
  p
}

#' make volcano plot from differential expression results
plot_volcano <- function(de, label = NULL) {
  p <- ggplot(de, aes(x = log2FoldChange, y = -log10(padj), color = -log10(padj))) +
    geom_vline(aes(xintercept = 0)) +
    geom_point(stroke = 0) +
    scale_color_gradient(low = "#f6edbd", high = "#008080") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    xlab(expression("Fold-change (log"[2]*")")) +
    ylab(expression(italic("P")*" value (-log"[10]*")")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text= element_text(size = 7, color = "black"),
          axis.title = element_text(size = 8))
  if (!is.null(label)) {
    p <- p +
      geom_text_repel(data = filter(de, gene_name %in% label),
                      aes(label = gene_name), color = "gray10", size = 2)
  }
  p
}

#' plot Sun1 cells on top of a given species UMAP projection
plot_umap_combined <- function(species = "mouse") {
  if (species == "mouse") reduction = "umap" else reduction = "ortho.umap"
  object <- get(species)
  df <- bind_rows(
    data.frame(
      "UMAP1" = object[[reduction]]@cell.embeddings[, 1],
      "UMAP2" = object[[reduction]]@cell.embeddings[, 2]
    ),
    data.frame(
      "UMAP1" = hypo[[paste0(species, ".umap")]]@cell.embeddings[, 1],
      "UMAP2" = hypo[[paste0(species, ".umap")]]@cell.embeddings[, 2],
      "Cluster" = Idents(hypo)
    )
  )
  # plot
  centers <- df %>%
    filter(Cluster %in% Idents(merged)) %>%
    group_by(Cluster) %>%
    summarize("UMAP1" = median(UMAP1), "UMAP2" = median(UMAP2))
  ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Cluster), size = 0.6, stroke = 0) +
    ggrepel::geom_text_repel(data = centers, aes(label = Cluster), size = 2) +
    scale_color_manual(values = colors, na.value = "gray90") +
    theme_void() +
    theme(legend.position = "none")
}

#' load saved cluster colors from a given species
load_colors <- function(species) {
  file <- paste0("results/conservation/", species, "_neuron_cluster_colors.csv")
  colors <- read.csv(file)
  colors <- tibble::deframe(colors)
  names(colors) <- levels(Idents(get(species)))
  colors
}

#' plot top n genes from a given species markers
plot_species_heatmap <- function(species, n = 10) {
  markers <- read.csv(paste0("results/conservation/", species, "_neuron_markers.csv"))
  genes <- markers %>% filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    group_by(gene) %>%
    slice(1) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    slice(1:10) %>%
    pull(gene)
  get(species) %>%
    ScaleData(features = genes, split.by = "Sample") %>%
    plot_heatmap(genes = genes, flipped = TRUE,
                 bars = TRUE, colors = get(paste0(species, "_colors"))) +
    theme(legend.position = "bottom")
}

get_species_gene_average <- function(genes, clusters = NULL) {
  # set up data.frame
  if (!all(genes %in% rownames(merged[["RNA"]]@scale.data))) {
    stop("Cannot find all genes in scale.data or meta.data")
  }
  if (length(genes) == 1) fn <- mean else fn <- Matrix::rowMeans
  if (is.null(clusters)) clusters <- levels(Idents(merged))
  species <- unique(merged$Species)
  df <- lapply(species, function(x) {
    sapply(clusters, function(y) {
      cells <- merged$Species == x & Idents(merged) == y
      fn(merged[["RNA"]]@scale.data[genes, cells])
    })
  })
  df <- lapply(df, as.data.frame)
  if (length(genes) == 1) {
    df <- lapply(df, setNames, "Z")
    df <- lapply(df, tibble::rownames_to_column, "Cluster")
    df <- lapply(df, function(x) mutate(x, "Gene" = genes))
  } else {
    df <- lapply(df, tibble::rownames_to_column, "Gene")
    df <- lapply(df, tidyr::pivot_longer, cols = -Gene, names_to = "Cluster", values_to = "Z")
  }
  df <- lapply(1:length(df), function(x) mutate(df[[x]], "Species" = species[x]))
  df <- bind_rows(df)
  df <- mutate(df, Gene = factor(Gene, levels = genes))
  df <- mutate(df, Cluster = factor(Cluster, levels = clusters))
  df <- mutate(df, Species = factor(Species, levels = c("Mouse", "Rat", "Macaque")))
  df
}

get_species_feature_average <- function(feature, clusters = NULL) {
  # set up data.frame
  if (!feature %in% colnames(merged@meta.data)) {
    stop("Feature is not in the meta.data slot")
  }
  if (is.null(clusters)) clusters <- levels(Idents(merged))
  species <- unique(merged$Species)
  df <- lapply(species, function(x) {
    sapply(clusters, function(y) {
      cells <- merged$Species == x & Idents(merged) == y
      mean(merged@meta.data[cells, feature])
    })
  })
  df <- lapply(df, as.data.frame)
  df <- lapply(df, setNames, "Expr")
  df <- lapply(df, tibble::rownames_to_column, "Cluster")
  df <- lapply(df, function(x) mutate(x, "Feature" = feature))
  df <- lapply(1:length(df), function(x) mutate(df[[x]], "Species" = species[x]))
  df <- bind_rows(df)
  df <- mutate(df, Cluster = factor(Cluster, levels = clusters))
  df <- mutate(df, Species = factor(Species, levels = c("Mouse", "Rat", "Macaque")))
  df
}


#' plot genes
plot_species_genes <- function(genes, clipped = FALSE, ...) {
  df <- get_species_gene_average(genes, ...)
  df <- mutate(df, Z = ifelse(Z > 1, 1, Z))
  df <- df %>% group_by(Gene, Species) %>% mutate(Z = round(Z/max(Z) * 100, 0))
  if (clipped) df <- mutate(df, Z = ifelse(Z > 0, 100, 0))
  palette <- c(
    colorRampPalette(c("#f6edbd", "#672044"))(101),
    colorRampPalette(c("#f6edbd", "#123f5a"))(101),
    colorRampPalette(c("#f6edbd", "#3d5941"))(101)
  )
  names(palette) <- paste(
    rep(c("Mouse", "Rat", "Macaque"), each = 101),
    rep(seq(0, 100), times = 3)
  )
  df <- mutate(df, Color = paste(Species, Z))
  df <- mutate(df, X = as.numeric(Gene) + as.numeric(Species)/3+0.33)
  ggplot(df, aes(x = X, y = Cluster, color = Color)) +
    geom_point() +
    scale_color_manual(values = palette, na.value = "#f6edbd") +
    theme_void() +
    theme(legend.position = "none",
          axis.text.y = element_text(hjust = 1, size = 8),
          axis.text.x = element_text(hjust = 1, vjust = 0.5, face = "italic",
                                     angle = 90, size = 8)) +
    scale_y_discrete(labels = levels(df$Cluster)) +
    scale_x_continuous(breaks = 2:(length(genes)+1), labels = genes)
}

#' plot an individual species in the merged UMAP space
plot_merged <- function(species) {
  sun1_cells <- merged$Sample == "Mouse_Sample_1776-AT-1"
  if (species == "Sun1") {
    cells <- sun1_cells
    colors <- colors
    clusters <- merged$Cluster[cells]
    species <- "LepRb-Sun1"
  } else {
    cells <- Cells(merged)[merged$Species == species & !sun1_cells]
    colors <- merged_colors
    clusters <- Idents(merged)[cells]
  }
  df <- bind_cols(
    as.data.frame(merged[["umap"]]@cell.embeddings[cells, ]),
    data.frame("Cluster" = clusters)
  )
  labels <- df %>% group_by(Cluster) %>%
    summarize("UMAP_1" = median(UMAP_1), "UMAP_2" = median(UMAP_2)) %>%
    filter(str_detect(Cluster, "^[A-Z]"))
  ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Cluster), stroke = 0, size = 0.8) +
    geom_text(data = labels, aes(label = Cluster), size = 2) +
    ggtitle(species) +
    theme_void() +
    scale_color_manual(values = colors, guide = "none") +
    theme(plot.title = element_text(hjust = 0.5, size = 8))
}

# get correlation between Sun1 and a given species
get_correlation <- function(species) {
  sun1_cells <- merged$Sample == "Mouse_Sample_1776-AT-1"
  cells <- Cells(merged)[merged$Species == species & !sun1_cells]

  df <- expand.grid(levels(Idents(merged)), unique(merged$Cluster[sun1_cells]))
  df$R <- apply(df, 1, function(x) {
    cells <- intersect(cells, Cells(merged)[Idents(merged) == x["Var1"]])
    sun1_cells <- intersect(Cells(merged)[sun1_cells], Cells(merged)[merged$Cluster == x["Var2"]])
    cor(
      Matrix::rowMeans(merged[["RNA"]]@scale.data[, cells]),
      Matrix::rowMeans(merged[["RNA"]]@scale.data[, sun1_cells])
    )
  })
  df$Var2 <- factor(df$Var2, levels = levels(Idents(hypo)))
  cluster_order <- df %>% group_by(Var1) %>%
    arrange(desc(R)) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(Var2)) %>%
    pull(Var1)
  df$Var1 <- factor(df$Var1, levels = cluster_order)
  df$Species <- species
  df
}

# plot correlation
plot_correlation <- function(species = c("Mouse", "Rat", "Macaque")) {
  df <- lapply(species, get_correlation)
  # get cluster levels by mean level from each species
  cluster_levels <- lapply(df, function(x) levels(x$Var1))
  unique_levels <- unlist(cluster_levels) %>% unique()
  ranks <- sapply(unique_levels, function(x) {
    median(sapply(cluster_levels, function(y) which(y == x)))
  })
  ranks <- sort(ranks)
  df <- bind_rows(df)
  df$Species <- factor(df$Species, levels = species)
  df$Var1 <- factor(df$Var1, levels = names(ranks))
  # plot
  ggplot(df, aes(x = Var2, y = Var1, fill = R)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080", name = "r") +
    theme_void() +
    facet_wrap(~ Species, ncol = 1) +
    theme(strip.text = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          legend.key.height = unit(0.02, "in"),
          legend.key.width = unit(0.1, "in"))
}

#' get coexpression frequency of 2 genes
get_coexpression <- function(object, gene1, gene2) {
  table(
    object[["RNA"]]@counts[gene1, ] > 0 & object[["RNA"]]@counts[gene2, ] > 0,
    Idents(object)
  ) %>%
    as.data.frame() %>%
    filter(Var1 == TRUE) %>%
    select(Var2, Freq)
}
