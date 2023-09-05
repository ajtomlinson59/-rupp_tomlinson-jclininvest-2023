# - Analysis of Lepr-Sun1 snRNA-seq data --------------------------------------
library(stringr)
library(Seurat)
library(dplyr)
source("scripts/snrna-seq_functions.R")
if (!dir.exists("results/sun1")) dir.create("results/sun1")

# read in all 10X data
dir <- list.dirs(path = "data/sun1")
dir <- dir[str_detect(dir, "filtered")]
mtx <- Read10X(dir)

# - Preprocessing -------------------------------------------------------------
# identify doublets using scrublet
reticulate::source_python("scripts/scrublet.py")
doublet_scores <- score_doublets(mtx)

# keep genes expressed in at least 5 cells
expressed <- apply(mtx, 1, function(x) sum(x > 0) > 4)
mtx <- mtx[expressed, ]

# create Seurat object
srt <- CreateSeuratObject(mtx)
rm(mtx)
srt$Doublet <- doublet_scores
srt$Sample <- "Sample_1776-AT-1"

# remove low expressing cells
srt <- subset(srt, nFeature_RNA > 600)

# - Clustering ----------------------------------------------------------------
srt <- normalize_data(srt)
srt <- reduce_dims(srt)

# predict cell types from published data by running CCA
hypo <- readRDS("~/Projects/published_data/GSE87544/GSE87544.rds")
hypo$Celltype <- ifelse(
  str_detect(Idents(hypo), "^[G,H,z]"), "Neuron", as.character(Idents(hypo))
)
hypo$Celltype <- str_replace_all(hypo$Celltype, "[0-9]$", "")
anchors <- FindTransferAnchors(hypo, srt, reduction = "cca")
labels <- TransferData(anchors, hypo$Celltype, weight.reduction = "cca")

# smooth out predictions by 10-NN
srt <- FindNeighbors(srt, dims = 1:srt[["pca"]]@misc$sig_pcs)
celltypes <- sapply(Cells(srt), function(x) {
  celltype <- labels[x, "predicted.id"]
  neighbors <- find_neighbors(srt, x, k = 10)
  scores <- colMeans(labels[neighbors, 2:(ncol(labels)-1)])
  scores <- sort(scores, decreasing = TRUE)
  str_extract(names(scores)[1], "(?<=score\\.).+")
})

# cluster to call celltypes
srt <- FindClusters(srt, resolution = 2)

# clusters with > 90% of one cell type are legit, rest are doublets
cluster_celltypes <- table(Idents(srt), celltypes)
cluster_celltypes <- apply(cluster_celltypes, 1, function(x) {
  if (max(x)/sum(x) >= 0.9) {
    names(x)[x == max(x)]
  } else {
    "Doublets"
  }
})
srt$Celltype <- recode(Idents(srt), !!!cluster_celltypes)

# save all cells
saveRDS(srt, "results/sun1/all_cells.rds")

# - Neurons -------------------------------------------------------------------
srt <- subset(srt, Celltype == "Neuron")
srt <- reduce_dims(srt)

# removing non-hypothalamic neurons
srt$Region <- case_when(
  Idents(srt) %in% c(3, 6, 11, 14, 17, 22, 31) ~ "TH",
  Idents(srt) == 23 ~ "VTA",
  Idents(srt) == 21 ~ "STR",
  TRUE ~ "HY"
)

# save all neurons
saveRDS(srt, "results/sun1/neurons.rds")

# keep just hypothalamus
srt <- subset(srt, Region == "HY")
srt <- reduce_dims(srt)
srt <- FindNeighbors(srt, dims = 1:srt[["pca"]]@misc$sig_pcs)
srt <- cluster(srt, return_object = TRUE)

# get cluster orders
cluster_order <- order_clusters(srt)
Idents(srt) <- factor(Idents(srt), levels = cluster_order, labels = 1:length(cluster_order))

# get markers for naming clusters
markers <- FindAllMarkers(srt, only.pos = TRUE)

# naming clusters
clusters <- c(
  `1` = "Nr5a1",
  `2` = "Irx5",
  `3` = "GLU1",
  `4` = "GABA",
  `5` = "Nts",
  `6` = "Pomc",
  `7` = "GLU2",
  `8` = "GLU3",
  `9` = "Foxb1",
  `10` = "KNDy",
  `11` = "Ghrh",
  `12` = "Opn5",
  `13` = "Eomes",
  `14` = "Pvalb",
  `15` = "Tbx19",
  `16` = "Prlh",
  `17` = "Agrp",
  `18` = "Glp1r"
)
Idents(srt) <- recode(Idents(srt), !!!clusters)
markers$cluster <- recode(markers$cluster, !!!clusters)
srt$Cluster <- as.character(Idents(srt))

# order clusters by Lepr expression
lepr <- sapply(split(srt[["RNA"]]@data["Lepr", ], Idents(srt)), mean)
cluster_order <- sort(lepr, decreasing = TRUE) %>% names()
Idents(srt) <- factor(Idents(srt), levels = cluster_order)
markers$cluster <- factor(markers$cluster, levels = cluster_order)

# save files
saveRDS(srt, "results/sun1/hypo.rds")
write.csv(markers, "results/sun1/markers_hypo.csv", row.names = FALSE, na = "")

# - TRAP enrichment -----------------------------------------------------------
# load in TRAP data
enrichment <- read.csv("results/GSE162603/enrichment.csv")
enriched <- filter(enrichment, padj < 0.05 & log2FoldChange > 0) %>% pull(gene_name)
enriched <- enriched[enriched %in% rownames(srt)]
de <- read.csv("results/GSE162603/de.csv") %>% arrange(padj)
de <- filter(de, padj < 0.05) %>% pull(gene_name)
de <- de[de %in% rownames(srt)]

markers <- FindAllMarkers(srt, features = enriched, logfc.threshold = 0, min.pct = 0)
write.csv(markers, "results/sun1/trap_enriched_markers.csv", row.names = FALSE, na = "")

# get PCA loading of genes that are DE by leptin treatment
mtx <- ScaleData(srt, features = de) %>% GetAssayData(slot = "scale.data")
pca <- prcomp(t(mtx))

# reverse sign if Stat3 is negative
if (pca$rotation["Stat3", 1] > 0) score = pca$x[, 1] else score = -pca$x[, 1]
srt$Leptin <- score

# save
saveRDS(srt, "results/sun1/hypo.rds")
