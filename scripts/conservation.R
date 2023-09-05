# - Hypothalamus datasets from mouse, rat, and macaque ------------------------

library(stringr)
library(Seurat)
library(dplyr)
source("scripts/snrna-seq_functions.R")
if (!dir.exists("results/conservation")) dir.create("results/conservation")

# - Read in data --------------------------------------------------------------
# solo rat
rat_dir <- "data/mbh/NovaA-219/Sample_683-AA-1/outs/filtered_feature_bc_matrix"
rat_mouse_dir <- "data/mbh/NovaA-230/Sample_871-AA-3/outs/filtered_feature_bc_matrix"
nhp_dir <- "data/GSE172203"

# read in all data
rat <- Read10X(rat_dir)
rat <- list(rat, parse_species(rat_mouse_dir, "rat"))
names(rat) <- c("Sample_683-AA-1", "Sample_871-AA-3")
mouse <- list(parse_species(rat_mouse_dir, "mouse"))
names(mouse) <- "Sample_871-AA-3"

# read in NHP data from GSE172203
nhp <- Matrix::readMM(paste0(nhp_dir, "/GSE172203_matrix.mtx.gz"))
rownames(nhp) <- read.table(paste0(nhp_dir, "/GSE172203_features.tsv.gz")) %>% pull(X1)
colnames(nhp) <- read.table(paste0(nhp_dir, "/GSE172203_barcodes.tsv.gz"), skip = 1) %>%
  pull(V1)
nhp_metadata <- read.csv(paste0(nhp_dir, "/GSE172203_metadata.csv.gz"))

# split by sample
nhp <- lapply(c("_1", "_2"), function(x) nhp[, str_detect(colnames(nhp), x)])
names(nhp) <- unique(nhp_metadata$Sample)

# - Process all species -------------------------------------------------------
reticulate::source_python("scripts/scrublet.py")
process_data <- function(dataset) {
  # keep shared & expressed genes
  common <- lapply(dataset, rownames) %>% unlist() %>% table() %>% .[. == length(dataset)] %>% names()
  dataset <- lapply(dataset, function(x) x[common, ])
  # keep genes expressed in at least 5 cells in all samples
  expression <- sapply(dataset, function(x) apply(x, 1, function(y) sum(y > 0) > 4))
  expressed <- apply(expression, 1, all)
  dataset <- lapply(dataset, function(x) x[expressed, ])
  # get doublet scores from scrublet
  doublet_scores <- lapply(dataset, score_doublets)
  sample_names <- rep(names(dataset), times = sapply(dataset, ncol))
  # create seurat object
  dataset <- CreateSeuratObject(do.call(cbind, dataset))
  dataset$Sample <- sample_names
  dataset$Doublet <- unlist(doublet_scores)
  dataset <- subset(dataset, nFeature_RNA > 600)
  dataset <- normalize_data(dataset, splitby = "Sample")
  dataset <- reduce_dims(dataset)
  return(dataset)
}

mouse <- process_data(mouse)
rat <- process_data(rat)
nhp <- process_data(nhp)

# add species info
mouse$Species <- "Mouse"
rat$Species <- "Rat"
nhp$Species <- "Macaque"

# - Celltype identification ---------------------------------------------------
# predict cell types from published data by running CCA
hypo <- readRDS("~/Projects/published_data/GSE87544/GSE87544.rds")
hypo$Celltype <- ifelse(
  str_detect(Idents(hypo), "^[G,H,z]"), "Neuron", as.character(Idents(hypo))
)
hypo$Celltype <- str_replace_all(hypo$Celltype, "[0-9]$", "")

#' Cluster cells and predict celltypes
get_celltypes <- function(object) {
  object <- cluster(object)
  if (object$Species[1] == "Rat") {
    object_stash <- object
    object <- to_mouse(object, "rnorvegicus_gene_ensembl")
  } else if (object$Species[1] == "Macaque") {
    object_stash <- object
    object <- to_mouse(object, "mmulatta_gene_ensembl")
  }
  anchors <- FindTransferAnchors(hypo, object, reduction = "cca")
  labels <- TransferData(anchors, hypo$Celltype, weight.reduction = "cca")
  celltypes <- sapply(Cells(object), function(x) {
    celltype <- labels[x, "predicted.id"]
    neighbors <- find_neighbors(object, x)
    scores <- colMeans(labels[neighbors, 2:(ncol(labels)-1)])
    scores <- sort(scores, decreasing = TRUE)
    str_extract(names(scores)[1], "(?<=score\\.).+")
  })
  celltype_table <- table(Idents(object), celltypes)
  # if majority of cells in a cluster have doublet score > 0.3, call doublets
  doublet_table <- table(Idents(object), object$Doublet > 0.3)
  doublet_clusters <- apply(doublet_table, 1, function(x) x["TRUE"] > x["FALSE"])
  cluster_celltypes <- sapply(rownames(celltype_table), function(x) {
    if (doublet_clusters[x]) {
      "Doublets"
    } else {
      colnames(celltype_table)[celltype_table[x, ] == max(celltype_table[x, ])][1]
    }
  })
  if (object$Species[1] != "Mouse") object <- object_stash
  object$Celltype <- recode(Idents(object), !!!cluster_celltypes)
  return(object)
}

mouse <- get_celltypes(mouse)
rat <- get_celltypes(rat)
nhp <- get_celltypes(nhp)

rm(hypo)

# - Remove doublets and low depth cells ---------------------------------------
for (species in c("mouse", "rat", "nhp")) {
  object <- get(species)
  doublets <- Cells(object)[object$Celltype == "Doublets" | object$Doublet > 0.3]
  if (length(doublets) > 0) {
    write.csv(doublets, paste0("results/conservation/", species, "_doublets.csv"),
              row.names = FALSE)
  }
  object <- subset(object, cells = Cells(object)[!Cells(object) %in% doublets])
  # remove low depth clusters based on < 2 SD compared to rest of cells in celltype
  for (celltype in unique(object$Celltype)) {
    depth <- object$nFeature_RNA[object$Celltype == celltype]
    depth <- (depth - mean(depth))/sd(depth)
    z <- sapply(split(depth, Idents(object)[names(depth)]), median)
    z <- z[!is.na(z)]
    if (any(z < -2)) {
      remove <- Cells(object)[Idents(object) %in% names(z)[z < -2]]
      write.csv(remove, paste0("results/conservation/low-depth_", species, "_",
                               str_to_lower(celltype), "s.csv"),
                row.names = FALSE)
      object <- subset(object, cells = Cells(object)[!Cells(object) %in% remove])
    }
  }
  assign(species, object)
}

# - Save each species ---------------------------------------------------------
# reprocess each
mouse <- reduce_dims(mouse)
rat <- reduce_dims(rat)
nhp <- reduce_dims(nhp)

# save all cells
saveRDS(mouse, "results/conservation/mouse.rds")
saveRDS(rat, "results/conservation/rat.rds")
saveRDS(nhp, "results/conservation/nhp.rds")

# - Neurons -------------------------------------------------------------------
mouse <- subset(mouse, Celltype == "Neuron") %>% reduce_dims() %>% cluster()
rat <- subset(rat, Celltype == "Neuron") %>% reduce_dims() %>% cluster()
nhp <- subset(nhp, Celltype == "Neuron") %>% reduce_dims() %>% cluster()

# order clusters
get_cluster_order <- function(object) {
  cluster_order <- order_clusters(object)
  Idents(object) <- factor(Idents(object), levels = cluster_order, labels = 1:length(cluster_order))
  object$Cluster <- Idents(object)
  object
}
mouse <- get_cluster_order(mouse)
rat <- get_cluster_order(rat)
nhp <- get_cluster_order(nhp)

# find markers
for (species in c("mouse", "rat", "nhp")) {
  get(species) %>%
    FindAllMarkers(only.pos = TRUE) %>%
    write.csv(paste0("results/conservation/", species, "_neuron_markers.csv"),
              row.names = FALSE)
}

# - Sun1 dataset projection ---------------------------------------------------
sun1 <- readRDS("results/sun1/hypo.rds")

# turn rat and nhp genes into mouse orthologs & embed all with UMAP model
process_species <- function(species) {
  if (species == "rat") dataset = "rnorvegicus_gene_ensembl"
  if (species == "nhp") dataset = "mmulatta_gene_ensembl"
  object <- to_mouse(get(species), dataset = dataset)
  object <- reduce_dims(object, return.model = TRUE)
  return(object)
}
mouse <- RunUMAP(mouse, dims = 1:mouse[["pca"]]@misc$sig_pcs, return.model = TRUE)
ortho_rat <- process_species("rat")
ortho_nhp <- process_species("nhp")

# add to original object
rat[["ortho.umap"]] <- ortho_rat[["umap"]]
nhp[["ortho.umap"]] <- ortho_nhp[["umap"]]

# embed new data & add to original
map_query <- function(species) {
  if (species %in% c("rat", "nhp")) {
    object <- get(paste0("ortho_", species))
  } else {
    object <- get(species)
  }
  anchors <- FindTransferAnchors(
    reference = object,
    query = sun1,
    dims = 1:object[["pca"]]@misc$sig_pcs,
    reduction = "cca"
  )
  sun1 <- MapQuery(
    anchorset = anchors,
    reference = object,
    query = sun1,
    refdata = "Cluster",
    reference.reduction = "cca",
    reduction.model = "umap",
    projectumap.args = list(
      reduction.name = paste0(species, ".umap"),
      reduction.key = paste0(species, "UMAP_")
    )
  )
  # change column names to match species
  col_names <- c(
    paste0("predicted.", species, ".Cluster.score"),
    paste0("predicted.", species, ".Cluster")
  )
  to_change <- which(str_detect(colnames(sun1@meta.data), "predicted.id"))
  if (length(to_change) < 2) stop("not enough metadata columns")
  colnames(sun1@meta.data)[to_change] <- col_names
  # return object
  return(sun1)
}

for (species in c("mouse", "rat", "nhp")) sun1 <- map_query(species)

# save results
saveRDS(mouse, "results/conservation/mouse_neurons.rds")
saveRDS(rat, "results/conservation/rat_neurons.rds")
saveRDS(nhp, "results/conservation/nhp_neurons.rds")
saveRDS(sun1, "results/sun1/hypo.rds")

# - Integration ---------------------------------------------------------------
# integrating all datasets using Harmony
sun1$Species <- "Mouse"
srt <- merge(sun1, list(mouse, ortho_rat, ortho_nhp))

# harmonize
srt$Sample <- paste(srt$Species, srt$Sample, sep = "_")
srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt, split.by = "Sample")
srt <- RunPCA(srt, npcs = 100)
srt[["pca"]]@misc$sig_pcs <- knee_test(srt)
srt <- harmonize(srt)
srt <- RunUMAP(srt, reduction = "harmony", dims = 1:srt[["pca"]]@misc$sig_pcs)
srt <- FindNeighbors(srt, reduction = "harmony", dims = 1:srt[["pca"]]@misc$sig_pcs)
srt <- FindClusters(srt)

cluster_order <- order_clusters(srt)
Idents(srt) <- factor(Idents(srt), levels = cluster_order, labels = 1:length(cluster_order))

saveRDS(srt, "results/conservation/integrated_neurons.rds")


# - Find Sun1 cluster correlates ----------------------------------------------
# define cells by their Sun1 neighbors (if none, )
sun1_samples <- srt$Sample == "Mouse_Sample_1776-AT-1"
sun1_table <- table(Idents(srt)[sun1_samples], srt$Cluster[sun1_samples])

# define clusters where 80% of Sun1 cells ended up in a single cluster as Sun1
matches <- apply(sun1_table, 2, function(x) max(x) > sum(x)*0.8)

# split clusters that have more than 1 Sun1 cluster defined to them
best <- apply(sun1_table[, matches], 2, function(x) which(x == max(x)))

# define clustering which maximizes Sun1 cohesiveness
get_scores <- function(idents) {
  cluster_table <- table(idents, obj$Cluster)
  bests <- apply(cluster_table, 2, function(x) names(x)[x == max(x)])
  if (any(duplicated(bests))) {
    return(NA)
  } else {
    sapply(1:ncol(cluster_table), function(x) {
      cluster_table[bests[x], x]/sum(cluster_table[, x])
    }) %>%
    mean()
  }
}

while (any(duplicated(best))) {
  to_split <- best[duplicated(best)] %>% unique()
  sig_pcs <- srt[["pca"]]@misc$sig_pcs
  resolutions <- seq(0.1, 1, by = 0.1)
  srt$New <- as.character(Idents(srt))
  for (cluster in to_split) {
    obj <- subset(srt, idents = cluster)
    sun1_options <- names(best)[best == cluster]
    obj$Cluster <- ifelse(obj$Cluster %in% sun1_options, obj$Cluster, NA)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:sig_pcs)
    clustering <- lapply(resolutions, function(x) {
      FindClusters(obj, resolution = x)@active.ident
    })
    scores <- sapply(clustering, get_scores)
    Idents(obj) <- clustering[scores == max(scores, na.rm = TRUE)]
    srt$New[Cells(obj)] <- paste(cluster, Idents(obj), sep = "_")
  }
  Idents(srt) <- factor(srt$New)
  # name clusters based on Sun1
  sun1_samples <- srt$Sample == "Mouse_Sample_1776-AT-1"
  sun1_table <- table(Idents(srt)[sun1_samples], srt$Cluster[sun1_samples])
  matches <- apply(sun1_table, 2, function(x) max(x) > sum(x)*0.8)
  best <- apply(sun1_table[, matches], 2, function(x) rownames(sun1_table)[x == max(x)])
}

# assign names
new_names <- names(best)
names(new_names) <- best
Idents(srt) <- recode(Idents(srt), !!!new_names)

# add names back to original Seurat objects (call `Conserved`)
mouse$Conserved <- as.character(Idents(srt)[srt$Sample == "Mouse_Sample_871-AA-3"])
rat$Conserved <- as.character(Idents(srt)[srt$Species == "Rat"])
nhp$Conserved <- as.character(Idents(srt)[srt$Species == "Macaque"])

# save objects
saveRDS(mouse, "results/conservation/mouse_neurons.rds")
saveRDS(rat, "results/conservation/rat_neurons.rds")
saveRDS(nhp, "results/conservation/nhp_neurons.rds")

# - Markers -------------------------------------------------------------------
markers <- pbapply::pblapply(levels(Idents(srt)), function(x) {
  FindConservedMarkers(srt, x, grouping.var = "Species", only.pos = TRUE) %>%
    mutate("cluster" = x)
})
markers <- bind_rows(markers)
markers$p_val_adj <- apply(markers, 1, function(x) {
  p_vals <- c(x["Macaque_p_val_adj"], x["Mouse_p_val_adj"], x["Rat_p_val_adj"])
  p_vals <- ifelse(p_vals < 10^-200, 10^-200, p_vals) %>% as.numeric()
  metap::logitp(p_vals)$p
})
markers <- tibble::rownames_to_column(markers, "gene") %>%
  mutate(gene = str_replace_all(gene, "\\.+[0-9]+$", ""))

write.csv(markers, "results/conservation/integrated_neuron_markers.csv", row.names = FALSE)

# - Leptin score --------------------------------------------------------------
de <- read.csv("results/GSE162603/de.csv")
de <- filter(de, gene_name %in% rownames(srt) & padj < 0.05) %>% pull(gene_name)

# get PCA loading of genes that are DE by leptin treatment
mtx <- ScaleData(srt, features = de, split.by = "Sample") %>%
  GetAssayData(slot = "scale.data")
pca <- prcomp(t(mtx))

# reverse sign if Stat3 is negative
if (pca$rotation["Stat3", 1] > 0) score = pca$x[, 1] else score = -pca$x[, 1]
srt$Leptin <- score

saveRDS(srt, "results/conservation/integrated_neurons.rds")

# - Lepr expression -----------------------------------------------------------
# get lepr expression across conserved populations relative to glia average
get_lepr <- function(species, gene = "Lepr") {
  object <- get(species)
  all_cells <- readRDS(paste0("results/conservation/", species, ".rds"))
  glia <- !all_cells$Celltype %in% c("Neuron", "Doublets")
  expr <- mean(all_cells[["RNA"]]@data[gene, glia])
  range <- max(all_cells[["RNA"]]@data[gene, ]) - min(all_cells[["RNA"]]@data[gene, ])
  result <- sapply(
    split(object[["RNA"]]@data[gene, ], object$Conserved),
    function(x) mean((x-expr)/range)
  )
  result[str_detect(names(result), "[A-Z]")]
}

mouse_lepr <- get_lepr("mouse")
rat_lepr <- get_lepr("rat")
nhp_lepr <- get_lepr("nhp", "LEPR")

lepr <- do.call(cbind, list(mouse_lepr, rat_lepr, nhp_lepr))
colnames(lepr) <- c("Mouse", "Rat", "Macaque")
clusters <- levels(Idents(sun1))[levels(Idents(sun1)) %in% rownames(lepr)]
lepr <- lepr[clusters, ]

write.csv(lepr, "results/conservation/lepr_expression.csv")
