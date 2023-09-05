# - Analyzing data hypothalamic Lepr TRAP-seq (Allison et al. 2018) -----------

library(dplyr)
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/GSE162603")) dir.create("results/GSE162603")

# - Preprocessing -------------------------------------------------------------
# read in data
counts <- read.csv("data/GSE162603/GSE162603_counts.csv.gz")
metadata <- read.csv("data/GSE162603/GSE162603_metadata.csv.gz")

# make counts into conventional data.frame
genes <- select(counts, starts_with("gene"))
counts <- select(counts, -gene_name) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_id")

# - Enrichment ----------------------------------------------------------------
# add pairing between bead and sup
metadata$Pair <- rep(LETTERS[1:(nrow(metadata)/2)], each = 2)
metadata$Cells <- factor(metadata$Cells, levels = c("Sup", "Bead"))

enrichment <- DESeq2::DESeqDataSetFromMatrix(
  counts[, metadata$Sample_ID], metadata, ~ Pair + Cells
) %>%
  DESeq2::DESeq() %>%
  DESeq2::results(name = "Cells_Bead_vs_Sup") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(genes, ., by = "gene_id")

# - Differential expression with Leptin ---------------------------------------
# only keep data from 10hr leptin treatment (and bead samples)
treatments <- c("10hrPBS", "10hrLeptin")
metadata <- filter(metadata, Cells == "Bead" & Treatment %in% treatments)
metadata$Treatment <- factor(metadata$Treatment, levels = treatments, labels = c("PBS", "Leptin"))
counts <- counts[, metadata$Sample_ID]

# run DESeq2
de <- DESeq2::DESeqDataSetFromMatrix(
  counts[, metadata$Sample_ID], metadata, ~ Mutant + Treatment
) %>%
  DESeq2::DESeq() %>%
  DESeq2::results(name = "Treatment_Leptin_vs_PBS") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(genes, ., by = "gene_id")

# - Save result ---------------------------------------------------------------
write.csv(enrichment, "results/GSE162603/enrichment.csv", row.names = FALSE, na = "")
write.csv(de, "results/GSE162603/de.csv", row.names = FALSE, na = "")
