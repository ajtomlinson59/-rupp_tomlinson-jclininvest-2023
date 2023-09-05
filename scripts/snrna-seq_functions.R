# - Various helper functions --------------------------------------------------
library(Seurat)
library(dplyr)
options(bitmapType = "cairo")

#' Split rat+mouse dual experiment by species
parse_species <- function(dir, species = "mouse") {
  mtx <- Read10X(dir)
  colnames(mtx) <- str_replace(colnames(mtx), "\\-[0-9]", "")
  genes <- str_extract(rownames(mtx), "(?<=rna\\_)(.+)$")
  calls <- read.csv(paste0(dir, "/../analysis/gem_classification.csv"))
  calls <- mutate(calls, barcode = str_replace(barcode, "\\-[0-9]", ""))
  if (species == "mouse") {
    keep_cells <- filter(calls, call == "mm10_3.0_premrna")$barcode
    keep_genes <- str_detect(rownames(mtx), "mm10")
  } else if (species == "rat") {
    keep_cells <- filter(calls, call == "Rnor_6.0_premrna")$barcode
    keep_genes <- str_detect(rownames(mtx), "Rnor")
  }
  rownames(mtx) <- genes
  mtx <- mtx[keep_genes, keep_cells]
  return(mtx)
}

#' Convert genes to orthologs (other species to mouse)
to_mouse <- function(object, dataset) {
  genes <- rownames(object)
  datasets <- list.files("data/ensembl")
  if (paste0(dataset, ".tsv") %in% datasets) {
    orthologs <- read.table(paste0("data/ensembl/", dataset, ".tsv"), header = TRUE, sep = "\t")
  } else {
    orthologs <- biomaRt::getLDS(
      attributes = "external_gene_name",
      filters = "external_gene_name",
      values = genes,
      mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
      attributesL = "external_gene_name",
      martL = biomaRt::useMart("ensembl", dataset = dataset)
    )
    write.table(orthologs, paste0("data/ensembl/", dataset, ".tsv"), sep = "\t", row.names = FALSE)
  }
  # POMC is not in the macaque orthologs file
  if (dataset == "mmulatta_gene_ensembl") {
    orthologs <- bind_rows(
      orthologs,
      data.frame("Gene.name" = "Pomc", "Gene.name.1" = "POMC")
    )
  }
  # Foxb1 is not in the rat orthologs file
  if (dataset == "rnorvegicus_gene_ensembl") {
    orthologs <- bind_rows(
      orthologs,
      data.frame("Gene.name" = "Foxb1", "Gene.name.1" = "Foxb1")
    )
  }
  # only keep 1:1 orthologs (unless there is a direct name match)
  matches <- str_to_lower(orthologs$Gene.name) == str_to_lower(orthologs$Gene.name.1)
  mouse_dups <- orthologs$`Gene.name`[duplicated(orthologs$`Gene.name`)]
  other_dups <- orthologs$`Gene.name.1`[duplicated(orthologs$`Gene.name.1`)]
  dups <- orthologs$Gene.name %in% mouse_dups | orthologs$Gene.name.1 %in% other_dups
  orthologs <- orthologs[!dups | matches, ]
  # convert relevant slots to ortholog genes
  for (x in c("counts", "data", "scale.data", "var.features")) {
    mtx <- slot(object[["RNA"]], x)
    if (is.null(dim(mtx))) {
      mtx <- mtx[mtx %in% orthologs$Gene.name.1]
      mtx <- orthologs$Gene.name[match(mtx, orthologs$Gene.name.1)]
    } else {
      mtx <- mtx[rownames(mtx) %in% orthologs$Gene.name.1, ]
      rownames(mtx) <- orthologs$Gene.name[match(rownames(mtx), orthologs$Gene.name.1)]
    }
    slot(object[["RNA"]], x) <- mtx
  }
  return(object)
}

#' Normalize data using `scran`
normalize_data <- function(object, splitby = NULL, return_object = TRUE,
                           verbose = FALSE) {
  # scran function
  run_scran <- function(mtx) {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list("counts" = mtx))
    if (verbose) message("Running dimension reduction ...")
    clusters <- scran::quickCluster(sce, min.size = 0)
    if (verbose) message("Computing sum factors ...")
    sce <- scran::computeSumFactors(sce, clusters = clusters)
    if (verbose) message("Log normalizing counts ...")
    sce <- scater::logNormCounts(sce)
    SingleCellExperiment::logcounts(sce)
  }
  # run for all data or each split
  if (is.null(splitby)) {
    mtx <- run_scran(object[["RNA"]]@counts)
  } else {
    mtx <- lapply(unique(object@meta.data[, splitby]), function(x) {
      message("Normalizing ", x, " ...")
      run_scran(object[["RNA"]]@counts[, object@meta.data[, splitby] == x])
    })
    mtx <- do.call(cbind, mtx)
    mtx <- mtx[, Cells(object)]
  }
  # return either mtx or object
  if (return_object) {
    object[["RNA"]]@data <- mtx
    return(object)
  } else {
    return(mtx)
  }
}

#' Choose PCs from elbow in scree plot
knee_test <- function(object, reduction = "pca") {
  n_pc = ncol(object[[reduction]]@cell.embeddings)
  total_var <- sum(object[[reduction]]@stdev^2)
  percent_var <- cumsum(object[[reduction]]@stdev^2)/total_var * n_pc
  diminishing <- which(percent_var - dplyr::lag(percent_var) < 1)
  return(min(diminishing) - 1)
}

#' Harmonize PCA in the object
harmonize <- function(object, metadata = "Sample") {
  harmony_embeddings <- harmony::HarmonyMatrix(
    object[["pca"]]@cell.embeddings[, 1:object[["pca"]]@misc$sig_pcs],
    object@meta.data[, metadata],
    do_pca = FALSE
  )
  object@reductions$harmony <- SeuratObject::CreateDimReducObject(
    embeddings = harmony_embeddings,
    key = "PC_"
  )
  return(object)
}

#' Run dimension reduction pipeline
reduce_dims <- function(object, reduction = "pca", ...) {
  object <- FindVariableFeatures(object)
  object <- ScaleData(object, split.by = "Sample")
  object <- RunPCA(object, npcs = 100)
  object[["pca"]]@misc$sig_pcs <- knee_test(object)
  if (reduction == "harmony") object <- harmonize(object)
  object <- RunUMAP(object, reduction = reduction,
                    dims = 1:object[["pca"]]@misc$sig_pcs, ...)
  return(object)
}

#' Get mean silhouette width for a given clustering
get_silhouette_width <- function(clusters, distances) {
  if (length(unique(clusters)) > 1) {
    sil <- cluster::silhouette(as.numeric(clusters), distances)
    return(mean(sil[, 3]))
  } else {
    NA
  }
}

#' Cluster cells in a Seurat object to maximize the mean silhouette width
cluster <- function(object, resolutions = seq(0.2, 2, by = 0.2),
                    reduction = "pca",
                    return_object = TRUE) {
  object <- FindVariableFeatures(object, verbose = FALSE)
  object <- ScaleData(object, split.by = "Sample", verbose = FALSE)
  object <- RunPCA(object, npcs = min(ncol(object)-1, 100), verbose = FALSE)
  sig_pcs <- knee_test(object, "pca")
  object[["pca"]]@misc$sig_pcs <- sig_pcs
  if (reduction == "harmony") object <- harmonize(object)
  object <- FindNeighbors(object, dims = 1:sig_pcs, reduction = reduction, verbose = FALSE)
  repeat {
    clusters <- lapply(resolutions, function(x) {
      obj <- try(FindClusters(object, resolution = x, verbose = FALSE), TRUE)
      if (class(obj) == "Seurat") Idents(obj)
    })
    clusters <- clusters[!sapply(clusters, is.null)]
    clusters <- do.call(cbind, clusters)
    if (length(clusters) == 0) break()
    if (ncol(clusters) == 1) break()
    if (ncol(object) > 10000) {
      cells <- sample(Cells(object), 10000)
    } else {
      cells <- Cells(object)
    }
    if (!is.matrix(object[[reduction]]@cell.embeddings)) return(NULL)
    distances <- cluster::daisy(object[[reduction]]@cell.embeddings[cells, 1:sig_pcs])
    widths <- apply(clusters[match(cells, Cells(object)), ], 2, function(x) {
      get_silhouette_width(x, distances)
    })
    if (all(is.na(widths))) {
      best_width <- 1
    } else {
      best_width <- which(widths == max(widths, na.rm = TRUE))[1]
    }
    if (best_width == length(widths)) {
      resolutions <- seq(max(resolutions), max(resolutions) + 2, by = 0.2)
    } else {
      break()
    }
  }
  # assign maximum silhouette width clustering
  if (length(clusters) == 0) {
    clusters <- rep("x", ncol(object))
  } else if (ncol(clusters) == 1) {
    clusters <- clusters[, 1]
  } else {
    clusters <- clusters[, best_width]
  }
  if (return_object) {
    Idents(object) <- clusters
    return(object)
  } else {
    return(clusters)
  }
}

#' Find neighbors from SNN graph
find_neighbors <- function(object, cell, k = 15) {
  snn <- object@graphs$RNA_snn[, cell]
  snn <- snn[snn > 0]
  snn <- sort(snn, decreasing = TRUE)
  names(snn)[1:k]
}


#' Run subclustering on 80% of cells within a celltype
subsample_cluster <- function(object, n, ...) {
  set.seed(n)
  cells <- sample(Cells(object), round(ncol(object)*0.8, 0))
  object <- subset(object, cells = cells)
  cluster(object, return_object = FALSE, reduction = reduction)
}


#' Convert subsample results to a matrix
convert_subsample_clusters <- function(subsample_results) {
  cells <- lapply(subsample_results, names) %>% unlist() %>% unique()
  results <- matrix(nrow = length(cells), ncol = length(subsample_results))
  rownames(results) <- cells
  for (i in 1:length(subsample_results)) {
    missing <- cells[!cells %in% names(subsample_results[[i]])]
    to_add <- rep(NA, length(missing)) %>% setNames(missing)
    iteration <- c(subsample_results[[i]], to_add)
    results[, i] <- iteration[cells]
  }
  return(results)
}

#' Get coclustering frequency of pairs of clusters (k and l)
#'
#' @param k cluster A
#' @param l cluster B
#' @param max_size maximum number of cells per cluster to include
#'
get_coclustering <- function(k, l, max_size = 500) {
  k_cells <- names(clusters)[clusters == k]
  l_cells <- names(clusters)[clusters == l]
  if (length(k_cells) > max_size) k_cells <- sample(k_cells, max_size)
  if (length(l_cells) > max_size) l_cells <- sample(l_cells, max_size)
  both <- unique(c(k_cells, l_cells))
  co_matrix <- matrix(0, nrow = length(both), ncol = length(both))
  pr_matrix <- matrix(0, nrow = length(both), ncol = length(both))
  rownames(co_matrix) <- colnames(co_matrix) <- both
  rownames(pr_matrix) <- colnames(pr_matrix) <- both
  result <- apply(subsample_clusters, 2, function(x) {
    x <- x[!is.na(x)]
    x = x[intersect(names(x), both)]
    pr_matrix[names(x), names(x)] <- pr_matrix[names(x), names(x)] + 1L
    for (i in unique(x)) {
      y <- names(x)[x == i]
      co_matrix[y, y] <- co_matrix[y, y] + 1L
    }
    mean(co_matrix[k_cells, l_cells]/pr_matrix[k_cells, l_cells], na.rm = TRUE)
  })
  mean(result)
}


#' Get average probability that a cell co-clustered with other cells of a given cluster
get_coclustering_cell <- function(cell) {
  probs <- sapply(unique(clusters), function(x) {
    x_cells <- names(clusters)[clusters == x]
    results <- apply(subsample_clusters, 2, function(y) {
      if (is.na(y[cell])) {
        NA
      } else {
        sum(y[cell] == y[x_cells], na.rm = TRUE)/sum(!is.na(y[x_cells]))
      }
    })
    mean(results, na.rm = TRUE)
  })
  unique(clusters)[probs == max(probs)]
}


#' Find closest other cluster in PCA space
find_neighboring_cluster <- function(object, cluster, reduction = "pca") {
  sig_pcs <- object[[reduction]]@misc$sig_pcs
  pca <- sapply(unique(Idents(object)), function(x) {
    colMeans(object[[reduction]]@cell.embeddings[Idents(object) == as.character(x), 1:sig_pcs])
  })
  colnames(pca) <- unique(Idents(object))
  distances <- dist(t(pca)) %>% as.matrix()
  distances <- distances[, as.character(cluster)] %>% .[. != 0]
  names(distances)[distances == min(distances)]
}

#' Get markers for either all cells or closest other cluster
get_markers <- function(object, cluster, other = NULL) {
  markers <- FindMarkers(
    object,
    ident.1 = cluster,
    ident.2 = other,
    logfc.threshold = 1,
    pseudocount.use = 0.01,
    min.pct = 0.2,
    only.pos = TRUE,
    verbose = FALSE
  )
  return(markers)
}

#' Calculate number of enriched genes
get_de_score <- function(markers) {
  if (nrow(markers) == 0) return(0)
  filter(markers, p_val_adj < 0.05) %>%
    filter(pct.1 > 2*pct.2) %>%
    filter(pct.2 < 0.2) %>%
    nrow()
}

#' Order clusters by mean scaled genes
order_clusters <- function(object, genes = NULL) {
  if (is.null(genes)) genes <- object[["RNA"]]@var.features
  if (!all(genes %in% rownames(object[["RNA"]]@scale.data))) {
    to_scale <- union(genes, rownames(object[["RNA"]]@scale.data))
    object <- ScaleData(object, features = to_scale, vars.to.regress = "Sample")
  }
  mtx <- sapply(unique(Idents(object)), function(x) {
    rowSums(object[["RNA"]]@scale.data[genes, Idents(object) == x])
  })
  colnames(mtx) <- unique(Idents(object))
  # cluster by hierarchical tree
  tree <- hclust(dist(t(mtx)))
  return(tree$labels[tree$order])
}

#' Plot labels from transfer anchors
plot_labels <- function(object, labels) {
  cells <- rownames(labels)
  df <- data.frame(
    "UMAP1" = object[["umap"]]@cell.embeddings[cells, 1],
    "UMAP2" = object[["umap"]]@cell.embeddings[cells, 2],
    "Cluster" = labels$predicted.id,
    "Score" = labels$prediction.score.max
  )
  centers <- df %>% group_by(Cluster) %>%
    summarize("UMAP1" = median(UMAP1), "UMAP2" = median(UMAP2))
  p1 <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Cluster), stroke = 0, size = 0.6) +
    geom_text(data = centers, aes(label = Cluster), size = 2) +
    theme_void() +
    theme(legend.position = "none")
  p2 <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Score), stroke = 0, size = 0.6) +
    scale_color_gradient(low = "gray90", high = "red4") +
    theme_void() +
    theme(legend.key.width = unit(0.02, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
  cowplot::plot_grid(p1, p2)
}

#' Get luminance for a given color
get_luminance <- function(color) {
  rgb <- col2rgb(color)
  lum <- 0.299*rgb["red", ] + 0.587*rgb["green", ] + 0.114*rgb["blue", ]
  unname(lum)
}

#' Get colors
get_colors <- function(n, lum_range = c(100, 200)) {
  colors <- grDevices::colors(distinct = TRUE)
  # get rid of excess grays
  colors <- colors[!str_detect(colors, "^gray[0-9]")]
  # keep colors within luminance range
  luminances <- sapply(colors, get_luminance)
  colors <- colors[luminances > lum_range[1] & luminances < lum_range[2]]
  while (length(colors) < n) {
    random_hex <- paste(sample(c(LETTERS[1:6], seq(0, 9)), 6, replace = TRUE), collapse = "")
    colors <- c(colors, paste0("#", random_hex))
  }
  # cluster to separate colors into distinct bins
  clusters <- kmeans(t(col2rgb(colors)), n)
  colors <- sapply(unique(clusters$cluster), function(x) {
    sample(colors[clusters$cluster == x], 1)
  })
  return(colors)
}

#' Order colors based on distance in Lab space
order_colors <- function(colors) {
  color_names <- colors
  colors <- sapply(colors, function(x) {
    if (str_detect(x, "^#")) {
      rgb <- colorspace::hex2RGB(x)@coords
    } else {
      col2rgb(x) %>% t()
    }
  }) %>% t()
  colnames(colors) <- c("R", "G", "B")
  # convert to Lab space
  colors <- convertColor(colors, from = "sRGB", to = "Lab")
  rownames(colors) <- color_names
  color_order <- hclust(dist(colors))
  color_order <- color_order$labels[color_order$order]
  return(color_order)
}

#' Plot samples in UMAP space
plot_samples <- function(object, colors = NULL, n_cells = 40000) {
  cells_plot <- sample(Cells(object), min(ncol(object), n_cells))
  p <- DimPlot(object, group.by = "Sample", cells = cells_plot) +
    theme_void() +
    theme(plot.title = element_blank())
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)
  return(p)
}
