
# Load Bioconductor packages used to import genomic files and draw coverage and
# annotation tracks. svglite is called explicitly when writing the SVG.
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

# The reference uses accession-like chromosome names rather than UCSC-style
# "chr" names, so disable Gviz's automatic UCSC name assumptions.
options(ucscChromosomeNames = FALSE)

# All paths are relative to the Coverage_Plots RStudio project directory.
genome_name <- "CpBGF"
gtf_file <- "../CpBGF_genome_V1.1.gtf"

# These CPM-normalized BigWigs are plotted in this order.
track_files <- c("../coord_bams_500_clean_unique/4_F-RBP_IP.bw", "../coord_bams_500_clean_unique/2_WT_Bunchgrass_IP.bw")
track_names <- c("FRBP_IP", "WT_IP")
track_colors <- c("#0000FF", "#0000FF")

# Genomic window and masking boundary for the cgd4_1910 locus.
chr <- "CP141121.1"
from <- 231120
to <- 233000
mask_from <- 232925
out_pdf <- "cgd4_1910.pdf"
out_svg <- sub("\\.pdf$", ".svg", out_pdf)

# Hide coverage at and after mask_from while retaining a continuous x-axis.
clip_right_edge <- function(gr, chr, from, to, mask_from) {
  if (length(gr) == 0L || is.null(mask_from) || mask_from > to) return(gr)

  overlap <- start(gr) < mask_from & end(gr) >= mask_from
  if (any(overlap)) {
    end(gr)[overlap] <- mask_from - 1L
    gr <- gr[width(gr) > 0]
  }

  gr <- gr[start(gr) < mask_from]
  zero_tail <- GRanges(chr, IRanges(max(mask_from, from), to), score = 0)
  c(gr, zero_tail)
}

# Import only the requested window, restrict to the target chromosome, and
# normalize empty inputs before applying the mask.
read_track <- function(file) {
  gr <- import(file, which = GRanges(chr, IRanges(from, to)))
  if (length(gr) && chr %in% seqlevels(gr)) {
    gr <- keepSeqlevels(gr, chr, pruning.mode = "coarse")
  }
  gr <- clip_right_edge(gr, chr, from, to, mask_from)
  if (length(gr) == 0L) {
    gr <- GRanges(chr, IRanges(from, to), score = 0)
  }
  gr
}

# Import the GTF and build the transcript annotation row below the coverage.
build_gene_track <- function() {
  gtf <- import(gtf_file)
  gtf <- keepSeqlevels(gtf, chr, pruning.mode = "coarse")

  # rtracklayer can expose the feature type under either metadata-column name.
  feature_col <- intersect(c("type", "feature"), names(mcols(gtf)))[1]
  if (!is.na(feature_col)) {
    gtf <- gtf[mcols(gtf)[[feature_col]] %in% c("transcript", "mRNA")]
  }

  track <- AnnotationTrack(
    range = gtf,
    genome = genome_name,
    chromosome = chr,
    name = "",
    feature = "transcript"
  )

  displayPars(track) <- list(
    shape = "arrow",
    col = "black",
    fill = "black",
    fill2 = "black",
    stacking = "dense",
    showFeatureId = FALSE,
    just.group = "above",
    arrowHeadWidth = 30,
    arrowHeadMaxWidth = 30,
    cex = 0.01
  )

  track
}

# Read both tracks and derive one shared y-axis from their scores.
track_ranges <- lapply(track_files, read_track)
all_scores <- unlist(lapply(track_ranges, function(gr) {
  score <- mcols(gr)$score
  score[is.na(score)] <- 0
  as.numeric(score)
}))

# A high quantile prevents one extreme bin from flattening the remainder.
ymax <- as.numeric(quantile(all_scores, probs = 0.995, na.rm = TRUE))
if (!is.finite(ymax) || ymax <= 0) ymax <- 1
ylim <- c(-max(0.05 * ymax, 0.1), ymax)

# Build one histogram DataTrack per sample with common scaling and styling.
coverage_tracks <- Map(function(gr, name, color) {
  track <- DataTrack(
    range = gr,
    data = "score",
    type = "histogram",
    name = name,
    genome = genome_name,
    chromosome = chr,
    ylim = ylim
  )
  displayPars(track) <- list(
    col = color,
    col.histogram = color,
    fill = color,
    lwd = 1.2,
    col.baseline = NA,
    lwd.baseline = 0
  )
  track
}, track_ranges, track_names, track_colors)

# Draw once through a function so PDF and SVG are identical.
draw_plot <- function() {
  plotTracks(
    c(coverage_tracks, list(build_gene_track())),
    from = from,
    to = to,
    main = sprintf("%s:%d-%d", chr, from, to),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    cex.main = 0.9,
    margin = 50
  )
}

# Write publication-friendly vector outputs.
pdf(out_pdf, width = 7, height = 5)
draw_plot()
dev.off()

svglite::svglite(out_svg, width = 7, height = 5)
draw_plot()
dev.off()
