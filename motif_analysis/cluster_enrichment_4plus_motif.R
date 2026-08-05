# This script is for analyzing whether cells annotated for highly expressing genes with 4 or more YBX2 motif belong to a specific cluster


library(ggrepel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)    
library(Seurat)
library(scales)
library(DT)

# This is the cparvum single cell seurat object
seurat_obj <- readRDS("all.rds")
print(seurat_obj)


# This is output from fimo, after coverting from t2t gene id into iowa2 gene id
gene_list_full <- readLines("fimo_iowa2_name.txt")
valid_genes_full <- intersect(gene_list_full, rownames(seurat_obj))

seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(valid_genes_full),
                             name = "Master_Score")
master_score_col <- "Master_Score1"
enriched_cells <- rownames(seurat_obj@meta.data)[
  seurat_obj@meta.data[[master_score_col]] > 0
]
seurat_obj$Is_Enriched <- rownames(seurat_obj@meta.data) %in% enriched_cells
print(paste("Enriched cells (master set):", length(enriched_cells)))

# These are derived from fimo output file
subgroup_files <- c("fimo_1motif_iowa2_name.txt",
                    "fimo_2motif_iowa2_name.txt",
                    "fimo_3motif_iowa2_name.txt",
                    "fimo_4plusmotif_iowa2_name.txt")

for (i in 1:4) {
  genes <- readLines(subgroup_files[i])
  valid_genes <- intersect(genes, rownames(seurat_obj))
  
  col_prefix <- paste0("Sub", i)
  seurat_obj <- AddModuleScore(seurat_obj,
                               features = list(valid_genes),
                               name = col_prefix)
  
  raw_score_col    <- paste0(col_prefix, "1")        
  masked_score_col <- paste0(col_prefix, "_Masked")  
  
  seurat_obj@meta.data[[masked_score_col]] <- ifelse(
    seurat_obj$Is_Enriched & seurat_obj@meta.data[[raw_score_col]] > 0,
    seurat_obj@meta.data[[raw_score_col]],
    NA
  )
}


scores <- seurat_obj@meta.data$Sub4_Masked
enrichment_threshold <- quantile(scores, 0.9, na.rm = TRUE) 
print(paste("Threshold for 'Enriched' set to:", round(enrichment_threshold, 4)))

seurat_obj$Is_4plus_Enriched <- !is.na(scores) & scores >= enrichment_threshold
print(table(seurat_obj$Is_4plus_Enriched))

seurat_obj$Cluster <- factor(seurat_obj$Cluster, levels = as.character(1:18))


p_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "Cluster",
                  pt.size = 1.5) +
  ggtitle("UMAP — Clusters") +
  xlim(-11, 11) + ylim(-14, 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text=element_text(size = 24), 
    legend.key.size = unit(1.5, 'cm'), 
    legend.text=element_text(size=24))+
  FontSize(x.title = 30, y.title = 30)+
  coord_fixed(ratio = 1)

print(p_umap)


meta <- seurat_obj@meta.data


cluster_colors <- hue_pal()(18)
names(cluster_colors) <- as.character(1:18)
cluster_summary <- meta %>%
  group_by(Cluster) %>%
  summarise(n_enriched = sum(Is_4plus_Enriched),
            n_total    = n(),
            .groups = "drop") %>%
  mutate(
    n_not_enriched = n_total - n_enriched,
    pct_enriched = n_enriched / n_total * 100,
    pct_not_enriched = 100 - pct_enriched
  ) %>%
  mutate(Cluster = factor(Cluster, levels = as.character(1:18))) %>%
  arrange(Cluster)
x_colors <- cluster_colors[as.character(1:18)]


count_long <- cluster_summary %>%
  select(Cluster, n_enriched, n_not_enriched) %>%
  pivot_longer(cols = c(n_enriched, n_not_enriched),
               names_to = "Status", values_to = "Count") %>%
  mutate(Status = factor(Status, levels = c("n_enriched", "n_not_enriched")))

p_bar_count <- ggplot(count_long, aes(x = Cluster, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.2, position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = c("n_enriched" = "blue",
                               "n_not_enriched" = "grey85"),
                    labels = c("4+ Motif Enriched", "Not Enriched"),
                    name = "") +
  geom_text(data = cluster_summary, 
            aes(x = Cluster, y = n_enriched, label = n_enriched),
            inherit.aes = FALSE, vjust = -0.5, color = "black", size = 4, fontface = "bold") +
  geom_text(data = cluster_summary, 
            aes(x = Cluster, y = n_total, label = n_total),
            inherit.aes = FALSE, vjust = -0.5, size = 4, color = "black") +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(cluster_summary$n_total) * 1.15)) +
  labs(title = "4+ Motif Enriched Cells per Cluster (Raw Count)",
       x = "Cluster", y = "Number of Cells") +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    plot.title  = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )




pct_long <- cluster_summary %>%
  select(Cluster, pct_enriched, pct_not_enriched) %>%
  pivot_longer(cols = c(pct_enriched, pct_not_enriched),
               names_to = "Status", values_to = "Pct") %>%
  mutate(Status = factor(Status, levels = c("pct_enriched", "pct_not_enriched")))

p_bar_pct <- ggplot(pct_long, aes(x = Cluster, y = Pct, fill = Status)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.2, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("pct_enriched" = "blue",
                               "pct_not_enriched" = "grey85"),
                    labels = c("4+ Motif Enriched", "Not Enriched"),
                    name = "") +
  geom_text(data = subset(pct_long, Status == "pct_enriched"),
            aes(label = sprintf("%.1f%%", Pct)),
            position = position_stack(vjust = 0.5, reverse = FALSE),
            color = "white", size = 4, fontface = "bold") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "4+ Motif Enriched Cells per Cluster (Proportion)",
       x = "Cluster", y = "Percentage of Cells (%)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    plot.title  = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )

# Printing out figure output
print(p_bar_count)
print(p_bar_pct)


# Fisher's Exact Test with BH Correction
fisher_results <- lapply(levels(meta$Cluster), function(cl) {
  in_cluster <- meta$Cluster == cl
  enriched   <- meta$Is_4plus_Enriched
  
  a <- sum(in_cluster & enriched)
  b <- sum(!in_cluster & enriched)
  c <- sum(in_cluster & !enriched)
  d <- sum(!in_cluster & !enriched)
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft_greater <- fisher.test(mat, alternative = "greater")$p.value 
  
  data.frame(
    Cluster      = cl,
    n_enriched   = a,
    n_cluster    = a + c,
    pct_enriched = round(a / (a + c) * 100, 2),
    p_value      = ft_greater
  )
})

fisher_df <- do.call(rbind, fisher_results)
fisher_df$p_adj_BH <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$signif <- cut(fisher_df$p_adj_BH,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", "ns"))

print(fisher_df)

# export to stats to html
datatable(fisher_df)

