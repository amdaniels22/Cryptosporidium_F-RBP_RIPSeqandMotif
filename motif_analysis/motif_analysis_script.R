# ==============================================================================
# Motif Analysis Script for scRNA-seq Data
# ==============================================================================
# This script analyzes Y-BOX binding motif enrichment in single-cell RNA-seq data.
#
# INPUT FILES REQUIRED (downloaded automatically if missing):
#   - all.rds (Seurat object with UMAP embeddings)
#     Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232438
#   - fimo_iowa2_name.txt                (~900 genes defining enriched cells)
#   - fimo_1motif_iowa2_name.txt         (Subset: 1 motif per CDS)
#   - fimo_2motif_iowa2_name.txt         (Subset: 2 motifs per CDS)
#   - fimo_3motif_iowa2_name.txt         (Subset: 3 motifs per CDS)
#   - fimo_4plusmotif_iowa2_name.txt     (Subset: 4+ motifs per CDS)
#   - piedata.tsv                        (For pie chart visualizations)
#
# OUTPUT:
#   - UMAP feature plots showing motif enrichment
#   - Density plots of module scores
#   - Pie charts of motif occurrence distribution
# ==============================================================================

# -----------------------------------------------------------------------------
# SECTION 0: AUTOMATED PACKAGE INSTALLATION
# -----------------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Install required packages
cat("=== Checking and Installing Required Packages ===\n")
install_if_missing("Seurat")
install_if_missing("ggplot2")
install_if_missing("ggrepel")
install_if_missing("dplyr")
install_if_missing("FactoMineR")
install_if_missing("factoextra")
install_if_missing("pheatmap")
install_if_missing("scales")
install_if_missing("corrplot")
install_if_missing("ggpubr")
cat("All packages ready!\n\n")

# -----------------------------------------------------------------------------
# SECTION 1: DATA DOWNLOAD AND PREPARATION
# -----------------------------------------------------------------------------

download_and_prepare_data <- function() {
  cat("=== Data Download and Preparation ===\n")

  # Download all.rds if missing
  if (!file.exists("all.rds")) {
    cat("all.rds not found. Downloading from NCBI GEO...\n")

    # Construct GEO download URL
    acc <- "GSE232438"
    file_name <- "GSE232438_seurat.object.crypto.atlas.all.rds.gz"
    url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file&file=%s",
                   acc, file_name)

    cat(sprintf("Downloading from: %s\n", url))

    # Download file
    download.file(url, destfile = file_name, mode = "wb")

    # Unzip
    cat("Unzipping file...\n")
    con <- gzfile(file_name, "rb")
    writeBin(readBin(con, raw(), n = 1e9), "all.rds")
    close(con)
    file.remove(file_name)

    # Find and rename the unzipped file
    unzipped_files <- list.files(pattern = "\\.rds$")
    if (length(unzipped_files) == 0) {
      stop("No .rds file found after unzipping!")
    }

    if (length(unzipped_files) > 1) {
      cat(sprintf("Warning: Multiple .rds files found. Using: %s\n", unzipped_files[1]))
    }

    source_file <- unzipped_files[1]
    if (!identical(source_file, "all.rds")) {
      cat(sprintf("Renaming %s to all.rds\n", source_file))
      file.rename(source_file, "all.rds")
    }

    # Clean up gz file
    file.remove(file_name)
    cat("Download and preparation complete!\n")
  } else {
    cat("all.rds already exists. Skipping download.\n")
  }

  # Check for other required files
  required_files <- c(
    "fimo_iowa2_name.txt",
    "fimo_1motif_iowa2_name.txt",
    "fimo_2motif_iowa2_name.txt",
    "fimo_3motif_iowa2_name.txt",
    "fimo_4plusmotif_iowa2_name.txt",
    "piedata.tsv"
  )

  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required files (please place in working directory):\n",
      paste("  -", missing_files, collapse = "\n")
    )
  }

  cat("All required files present!\n\n")
}

# Run data preparation
download_and_prepare_data()

# -----------------------------------------------------------------------------
# SECTION 2: LOAD DATA
# -----------------------------------------------------------------------------

cat("=== Loading Data ===\n")

# Load required libraries (already installed in Section 0)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(pheatmap)
  library(scales)
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(ggpubr)
})

# Load Seurat object
seurat_obj <- readRDS("all.rds")
cat(sprintf("Loaded Seurat object with %d cells and %d features\n",
            ncol(seurat_obj), nrow(seurat_obj)))

# Read gene lists
cat("Reading gene lists...\n")
gene_list_full <- readLines("fimo_iowa2_name.txt")
cat(sprintf("Full gene list: %d genes\n", length(gene_list_full)))

subgroup_files <- c(
  "fimo_1motif_iowa2_name.txt",
  "fimo_2motif_iowa2_name.txt",
  "fimo_3motif_iowa2_name.txt",
  "fimo_4plusmotif_iowa2_name.txt"
)

# -----------------------------------------------------------------------------
# SECTION 3: MASTER SCORE - DEFINE ENRICHED CELL SET
# -----------------------------------------------------------------------------

cat("\n=== Computing Master Score ===\n")

# Validate genes exist in Seurat object
valid_genes_full <- intersect(gene_list_full, rownames(seurat_obj))
cat(sprintf("Valid genes (in Seurat object): %d\n", length(valid_genes_full)))

# Compute module score for all genes
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(valid_genes_full),
  name = "Master_Score"
)

master_score_col <- "Master_Score1"
cat(sprintf("Created metadata column: %s\n", master_score_col))

# Define enriched cells: those with Master_Score > 0
enriched_cells <- rownames(seurat_obj@meta.data)[
  seurat_obj@meta.data[[master_score_col]] > 0
]
seurat_obj$Is_Enriched <- rownames(seurat_obj@meta.data) %in% enriched_cells

cat(sprintf("Enriched cells (master set): %d (%.1f%%)\n",
            length(enriched_cells),
            100 * length(enriched_cells) / ncol(seurat_obj)))

# -----------------------------------------------------------------------------
# SECTION 4: SUB-GROUP SCORES WITH MASKING
# -----------------------------------------------------------------------------

cat("\n=== Computing Sub-group Scores ===\n")

# Store sub-group info
subgroup_info <- data.frame(
  name = character(),
  n_genes = integer(),
  score_col = character(),
  stringsAsFactors = FALSE
)

for (i in 1:4) {
  genes <- readLines(subgroup_files[i])
  valid_genes <- intersect(genes, rownames(seurat_obj))

  col_prefix <- paste0("Sub", i)
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(valid_genes),
    name = col_prefix
  )

  raw_score_col <- paste0(col_prefix, "1")
  masked_score_col <- paste0(col_prefix, "_Masked")

  # Mask: only keep score if cell is in enriched set AND sub-group score > 0
  seurat_obj@meta.data[[masked_score_col]] <- ifelse(
    seurat_obj$Is_Enriched & seurat_obj@meta.data[[raw_score_col]] > 0,
    seurat_obj@meta.data[[raw_score_col]],
    NA
  )

  subgroup_info <- rbind(subgroup_info, data.frame(
    name = c("1 Motif", "2 Motifs", "3 Motifs", "4+ Motifs")[i],
    n_genes = length(valid_genes),
    score_col = masked_score_col,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("Subgroup %d: %d genes -> %s\n",
              i, length(valid_genes), masked_score_col))
}

print(subgroup_info)

# -----------------------------------------------------------------------------
# SECTION 5: UMAP FEATURE PLOTS - MASTER SCORE
# -----------------------------------------------------------------------------

cat("\n=== Generating UMAP Plots ===\n")

plot_umap <- function(seurat_obj, features, title, xlim = c(-11, 11),
                      ylim = c(-14, 8), pt.size = 2) {

  p <- FeaturePlot(
    object = seurat_obj,
    features = features,
    reduction = "umap",
    order = TRUE,
    pt.size = pt.size,
    cols = c("lightgrey", "grey80", "blue")
  ) +
    ggtitle(title) +
    xlim(xlim) + ylim(ylim) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text = element_text(size = 14),
      legend.key.size = unit(1, 'cm'),
      legend.text = element_text(size = 12)
    ) +
    coord_fixed(ratio = 1)

  return(p)
}

# Master Score
p_master <- plot_umap(
  seurat_obj,
  "Master_Score1",
  "UMAP: Y-BOX Binding Motif Gene Set Enrichment"
)
print(p_master)

# Sub-group plots with masked scores
for (i in 1:4) {
  masked_col <- paste0("Sub", i, "_Masked")
  title <- c("1 Motif", "2 Motifs", "3 Motifs", "4+ Motifs")[i]

  p <- plot_umap(
    seurat_obj,
    masked_col,
    paste("UMAP:", title, "Found per CDS"),
    pt.size = 1.5
  ) +
    scale_color_gradientn(
      colours = c("lightgrey", "grey80", "blue"),
      na.value = "lightgrey"
    )

  print(p)
}

# -----------------------------------------------------------------------------
# SECTION 6: MODULE SCORE DENSITY PLOT
# -----------------------------------------------------------------------------

cat("\n=== Computing Density Plot ===\n")

module_scores <- seurat_obj@meta.data$Master_Score1
d <- density(module_scores, bw = "SJ", n = 512)

# Find valley between peaks (cutoff point)
valley_region <- which(d$x > 0 & d$x < 0.1)
valley_idx <- valley_region[which.min(d$y[valley_region])]
valley_cutoff <- d$x[valley_idx]

cat(sprintf("Density valley cutoff: %.4f\n", valley_cutoff))

density_plot <- ggplot(data.frame(score = module_scores), aes(x = score)) +
  geom_density(fill = "lightblue", alpha = 0.5, size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = valley_cutoff, linetype = "dashed", color = "black",
             linewidth = 0.8) +
  annotate("text", x = valley_cutoff, y = max(d$y) * 0.9,
           label = sprintf("Enrichment Cutoff: %.4f", valley_cutoff),
           color = "black", hjust = -0.1, size = 5) +
  ggtitle("Master Score Density Distribution") +
  labs(x = "Module Score", y = "Density") +
  theme_bw()

print(density_plot)

# -----------------------------------------------------------------------------
# SECTION 7: PIE CHARTS - MOTIF OCCURRENCE DISTRIBUTION
# -----------------------------------------------------------------------------

cat("\n=== Generating Pie Charts ===\n")

# Read pie chart data
piedata <- read.table("piedata.tsv", header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
colnames(piedata) <- c("ID", "Value")

# Create categories
piedata$Category <- "Value > 3"
piedata$Category[piedata$Value == 1] <- "Value = 1"
piedata$Category[piedata$Value == 2] <- "Value = 2"
piedata$Category[piedata$Value == 3] <- "Value = 3"

# Count frequencies
count_data <- as.data.frame(table(piedata$Category))
colnames(count_data) <- c("Category", "Count")

# Ensure all categories exist in specified order
all_cats <- c("Value = 1", "Value = 2", "Value = 3", "Value > 3")
plot_data <- merge(data.frame(Category = all_cats, stringsAsFactors = FALSE),
                   count_data, all.x = TRUE)
plot_data$Count[is.na(plot_data$Count)] <- 0
plot_data <- plot_data[match(all_cats, plot_data$Category), ]

# Pie chart 1: Highlight Value = 1
colors_1 <- c("skyblue", "lightgrey", "lightgrey", "lightgrey")
pie(plot_data$Count, labels = plot_data$Category, col = colors_1,
    main = "Y-BOX Motifs per Gene = 1", border = "white")

# Pie chart 2: Highlight Value = 2
colors_2 <- c("lightgrey", "skyblue", "lightgrey", "lightgrey")
pie(plot_data$Count, labels = plot_data$Category, col = colors_2,
    main = "Y-BOX Motifs per Gene = 2", border = "white")

# Pie chart 3: Highlight Value = 3
colors_3 <- c("lightgrey", "lightgrey", "skyblue", "lightgrey")
pie(plot_data$Count, labels = plot_data$Category, col = colors_3,
    main = "Y-BOX Motifs per Gene = 3", border = "white")

# Pie chart 4: Highlight Value > 3
colors_4 <- c("lightgrey", "lightgrey", "lightgrey", "skyblue")
pie(plot_data$Count, labels = plot_data$Category, col = colors_4,
    main = "Y-BOX Motifs per Gene > 3", border = "white")

# Summary pie chart with labels
freq_table <- table(piedata$Category)
freq_table <- freq_table[all_cats]  # Ensure correct order
labels_with_counts <- paste0(plot_data$Category, " (n=", plot_data$Count, ")")

pie(freq_table, labels = labels_with_counts,
    main = "Motif Occurrence Distribution",
    col = grey.colors(length(freq_table), start = 0.9, end = 0.3),
    border = "white")

# -----------------------------------------------------------------------------
# SECTION 8: SUMMARY
# -----------------------------------------------------------------------------

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Enriched cells: %d (%.1f%%)\n",
            length(enriched_cells),
            100 * length(enriched_cells) / ncol(seurat_obj)))
cat("Output: UMAP plots and pie charts generated\n")
