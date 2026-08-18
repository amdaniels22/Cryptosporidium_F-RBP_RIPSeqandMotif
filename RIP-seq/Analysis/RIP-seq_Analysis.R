
############### ############### ############### ############### ###############
############### This R script uses salmon generated TPM quantification of 
############### RIP-seq data as input, performs analysis using that data
############### in order to identify putative F-RBP target transcripts
############### and produces Figure 5B, Suppl. Figure 5,
############### as well as Suppl. Table 2
############### ############### ############### ############### ###############

library(ggplot2)
library(mixtools)
library(tximport)
library(svglite)

############### ############### ############### ############### ###############
############### import salmon quant outputs (TPM estimates)
############### ############### ############### ############### ###############

salmon_dir <- "../salmon_quant"
tx2gene_file <- "tx2gene.csv"

# Folder names produced by the Salmon command in the top-level README.
sample_dirs <- c(
  BG_Input = "1_WT_Bunchgrass_Input.trim",
  BG_IP = "2_WT_Bunchgrass_IP.trim",
  FRBP_Input = "3_F-RBP_Input.trim",
  FRBP_IP = "4_F-RBP_IP.trim"
)

files <- file.path(salmon_dir, unname(sample_dirs), "quant.sf")
names(files) <- names(sample_dirs)

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop(
    "Missing Salmon quantification file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\nOpen Analysis.Rproj before running this script."
  )
}

if (!file.exists(tx2gene_file)) {
  stop("Missing transcript-to-gene table: ", tx2gene_file)
}

tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)
if (!identical(colnames(tx2gene)[1:2], c("transcript_id", "gene_id"))) {
  stop("tx2gene.csv must begin with columns transcript_id and gene_id.")
}
if (anyDuplicated(tx2gene$transcript_id)) {
  stop("tx2gene.csv contains duplicated transcript IDs.")
}

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"   # recommended
)

TPM <- as.data.frame(txi$abundance)
colnames(TPM) <- paste0(colnames(TPM), "_TPM")
colSums(TPM)
TPM$Name <- rownames(TPM)

Counts <- as.data.frame(txi$counts)
colSums(Counts)
Counts$Name <- rownames(Counts)

data <- merge(TPM, Counts, by = "Name")

# Keep a fixed, documented column order independent of filesystem ordering.
expected_columns <- c(
  "Name",
  "BG_Input_TPM", "BG_IP_TPM", "FRBP_Input_TPM", "FRBP_IP_TPM",
  "BG_Input", "BG_IP", "FRBP_Input", "FRBP_IP"
)
data <- data[, expected_columns]

### only include genes with a length-scaled estimated count greater than 5 in each sample

data_all <- data 

data <- data[data$BG_Input>5&data$BG_IP>5&data$FRBP_Input>5&data$FRBP_IP>5,]


############### ############### ############### ############### ###############
############### calculate different ratios and their logs
############### ############### ############### ############### ###############

data$BG_ratio<-(data$BG_IP_TPM/data$BG_Input_TPM)
data$BG_log_ratio<-log2(data$BG_IP_TPM/data$BG_Input_TPM)

data$FRBP_ratio<-(data$FRBP_IP_TPM/data$FRBP_Input_TPM)
data$FRBP_log_ratio<-log2(data$FRBP_IP_TPM/data$FRBP_Input_TPM)

data$enriched<-log2(data$FRBP_ratio/data$BG_ratio)

### get a feel for the distribution of some of those ratios by calculating their means

mean(data$BG_log_ratio)
mean(data$FRBP_log_ratio)
mean(data$enriched)

############### ############### ############### ############### ###############
############### mixtools for mixture modeling
############### ############### ############### ############### ###############

set.seed(123)  # for reproducibility
mix_model <- normalmixEM(data$enriched,
                         k = 2,
                         maxit = 100000,
                         arbvar = F)  # enforce equal variance

# Identify components by their fitted means so component-label switching cannot
# reverse the enriched/non-enriched interpretation.
background_component <- which.min(mix_model$mu)
enriched_component <- which.max(mix_model$mu)

data$posterior_background <- mix_model$posterior[, background_component]
data$posterior_enriched <- mix_model$posterior[, enriched_component]

# Compute Log Odds Ratio (LOD)
data$LOD <- log(data$posterior_enriched / data$posterior_background)

# Classify enriched genes
data$model_enriched <- data$LOD > 0

# How many non-enriched genes vs. enriched genes ? Expect mostly FALSE, smaller sub-population of F-RBP targets = TRUE
table(data$model_enriched)

############### ############### ############### ############### ###############
############### plotting of salmon TPM data and mixture modeling results
############### ############### ############### ############### ###############

dir.create("output", showWarnings = FALSE, recursive = TRUE)

### Suppl. Figure 5: histogram of "enriched" column and overlay line plots of the two mixture model distributions
hist(data$enriched, breaks = 50, probability = T, main = "Mixture Model", border = F, xlab = "log2((FRBP_IP_TPM/FRBP_Input_TPM)/(BG_IP_TPM/BG_Input_TPM))")
curve(mix_model$lambda[background_component] *
        dnorm(x, mix_model$mu[background_component], mix_model$sigma[background_component]),
      col = "black", add = TRUE, lwd = 4)
curve(mix_model$lambda[enriched_component] *
        dnorm(x, mix_model$mu[enriched_component], mix_model$sigma[enriched_component]),
      col = "#0000FF", add = TRUE, lwd = 4)


svglite::svglite("./output/mixture_model.svg", width = 6, height = 4)

hist(data$enriched, breaks = 50, probability = TRUE,
     main = "Mixture Model", border = FALSE,
     xlab = "log2((FRBP_IP_TPM/FRBP_Input_TPM)/(BG_IP_TPM/BG_Input_TPM))")

curve(mix_model$lambda[background_component] *
        dnorm(x, mix_model$mu[background_component], mix_model$sigma[background_component]),
      col = "black", add = TRUE, lwd = 3)

curve(mix_model$lambda[enriched_component] *
        dnorm(x, mix_model$mu[enriched_component], mix_model$sigma[enriched_component]),
      col = "#0000FF", add = TRUE, lwd = 3)

dev.off()

### Figure 5B: plot BG and FRBP IP/Input log ratios of each gene and highlight genes predicted to be F-RBP targets based on mixture model
# note: 1 data point was designated as an outlier and excluded from plot for improved readability

ggplot(data = data, aes(BG_log_ratio, FRBP_log_ratio, color = model_enriched&FRBP_log_ratio>0.5))+geom_point()+
  scale_color_manual(values = c("TRUE" = "#0000FF", "FALSE" = "gray")) +
  theme_classic() +  # white background
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # black border
    axis.line = element_blank(), # removes default axis lines
    legend.position = "none"  # removes legend
  ) +
  coord_fixed() +  # same scale aspect ratio for x and y
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1))  

p <- ggplot(data = data, aes(BG_log_ratio, FRBP_log_ratio, color = model_enriched&FRBP_log_ratio>0.5))+geom_point()+
  scale_color_manual(values = c("TRUE" = "#0000FF", "FALSE" = "gray")) +
  theme_classic() +  # white background
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # black border
    axis.line = element_blank(), # removes default axis lines
    legend.position = "none"  # removes legend
  ) +
  coord_fixed() +  # same scale aspect ratio for x and y
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1))  

ggsave("./output/BG_vs_FRBP_log_ratio.svg", plot = p, device = svglite::svglite, width = 3, height = 3)

############### ############### ############### ############### ###############
############### create new data frame with only F-RBP target transcripts
############### = basis for Suppl. Table 2
############### ############### ############### ############### ###############

hits<-data[data$model_enriched==T&data$FRBP_log_ratio>0.5,]
### Note: only transcripts that are predicted targets based on mixture modeling 
### AND that show sufficiently increased TP in F-RBP_IP vs. F-RBP input are considered 
### in final F-RBP target list

CpBGF <- read.csv(file = "CpBGF_gene_table.csv")
colnames(CpBGF)[7]<-"Name"

hits <- merge(hits, CpBGF, by = "Name")

gamete_markers <- read.csv(file = "Walzer_2024_sex_specific_genes.csv")
gamete_markers <- gamete_markers[,1:5]

hits$sex_specific <- gamete_markers$sex[match(hits$old_ID, gamete_markers$geneID)]
hits$sex_specific[ is.na(hits$sex_specific) ] <- "none"

hits <- hits[, !names(hits) %in% c("model_enriched", "gbkey", "Name_Conversion", "Notes")]

write.csv(hits, file = "./output/RIP-seq_hits.csv", row.names = FALSE)
