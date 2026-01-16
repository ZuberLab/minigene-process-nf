#!/usr/bin/env Rscript

################################################################################
# perform PCA on count data
# part of minigene screen pre-processing pipeline
#
# Florian Andersch
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2026/01/16
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(plotly)
})

args <- commandArgs(trailingOnly = TRUE)

# ============================
# User settings (edit these)
# ============================

input_file <- args[1]     # first two columns: id, group; remaining columns: samples
outdir     <- "pca"
sep        <- "\t"
header     <- TRUE

# Transform / PCA settings
prior_count     <- 1
topN            <- 30
topN_label <- 10   # label only the top N loading arrows (still draw topN arrows)
center_features <- TRUE
scale_features  <- FALSE

# Metadata (Excel). Set to NULL to disable.
metadata_xlsx   <- args[2]      # must contain column 'sample_id' matching count-table sample column names
metadata_sheet  <- 1
sample_id_col   <- "sample_id"
color_by        <- "condition"
shape_by        <- "cellline_name"

# Biplot arrow scaling (0.7–1.0 usually looks good)
arrow_scale_fraction <- 0.8

# ============================
# Helpers
# ============================

stop_if <- function(cond, msg) {
  if (cond) stop(msg, call. = FALSE)
}

# Median ratio normalization (DESeq-style) + log2
log_mrn <- function(counts, prior_count = 1) {
  stop_if(any(counts < 0, na.rm = TRUE), "Counts contain negative values.")
  stop_if(!is.matrix(counts), "counts must be a numeric matrix.")
  
  log_counts <- log(counts)
  log_counts[!is.finite(log_counts)] <- NA  # log(0) -> NA
  
  gm_log <- rowMeans(log_counts, na.rm = TRUE)
  gm <- exp(gm_log)
  
  keep <- is.finite(gm) & gm > 0
  stop_if(sum(keep) < 2, "Too few non-zero features to compute size factors.")
  
  ratios <- sweep(counts[keep, , drop = FALSE], 1, gm[keep], FUN = "/")
  
  size_factors <- apply(ratios, 2, function(x) {
    x <- x[is.finite(x) & x > 0]
    stop_if(length(x) == 0, "A sample has no finite positive ratios; cannot compute size factor.")
    median(x)
  })
  
  stop_if(any(!is.finite(size_factors)) || any(size_factors <= 0),
          "Invalid size factors computed.")
  
  norm_counts <- sweep(counts, 2, size_factors, FUN = "/")
  log2(norm_counts + prior_count)
}

# ============================
# Main
# ============================

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Reading input: ", input_file)
df <- read.table(
  input_file,
  header = header,
  sep = sep,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stop_if(ncol(df) < 4, "Expected at least 4 columns: id, group, and >=2 sample columns.")

id_col <- 1
group_col <- 2

feature_ids <- df[[id_col]]
stop_if(anyDuplicated(feature_ids) > 0, "Feature IDs in 'id' column contain duplicates. Make them unique first (or aggregate).")

counts_df <- df[, -(c(id_col, group_col)), drop = FALSE]

counts_mat <- as.matrix(counts_df)
suppressWarnings(storage.mode(counts_mat) <- "numeric")
stop_if(any(is.na(counts_mat)), "Counts matrix contains NA after numeric conversion. Check input formatting.")
rownames(counts_mat) <- feature_ids

message("Normalizing: median ratio normalization (MRN) + log2")
log_mat <- log_mrn(counts_mat, prior_count = prior_count)

message("Running PCA")
pca <- prcomp(t(log_mat), center = center_features, scale. = scale_features)

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
var_df <- data.frame(
  PC = paste0("PC", seq_along(var_expl)),
  variance_explained = var_expl,
  stringsAsFactors = FALSE
)

scores <- as.data.frame(pca$x)
scores$sample_id <- rownames(scores)

loadings <- as.data.frame(pca$rotation)
loadings$feature <- rownames(pca$rotation)

# Save loadings + variance explained now (independent of metadata)
write.table(loadings, file.path(outdir, "pca_loadings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(var_df,   file.path(outdir, "pca_variance_explained.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------
# Read metadata from Excel (optional) and merge by sample_id
# ----------------------------
if (!is.null(metadata_xlsx)) {
  message("Reading metadata (Excel): ", metadata_xlsx)
  meta <- readxl::read_excel(metadata_xlsx, sheet = metadata_sheet)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  
  stop_if(!(sample_id_col %in% colnames(meta)),
          paste0("Metadata Excel must contain a column named '", sample_id_col, "'."))
  
  colnames(meta)[colnames(meta) == sample_id_col] <- "sample_id"
  
  missing_meta <- setdiff(scores$sample_id, meta$sample_id)
  if (length(missing_meta) > 0) {
    warning("Some samples in counts are missing from metadata: ",
            paste(missing_meta, collapse = ", "))
  }
  
  scores <- merge(scores, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
}

# Save scores AFTER metadata merge (so TSV includes metadata columns)
write.table(scores, file.path(outdir, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", file.path(outdir, "pca_scores.tsv"))
message("Wrote: ", file.path(outdir, "pca_loadings.tsv"))
message("Wrote: ", file.path(outdir, "pca_variance_explained.tsv"))

# ----------------------------
# Plot: PCA PC1 vs PC2
# ----------------------------
pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

p_base <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_point(size = 3, alpha = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, title = "PCA: PC1 vs PC2") +
  theme_bw(base_size = 12)

if (!is.null(color_by)) {
  stop_if(is.null(metadata_xlsx), "color_by is set but metadata_xlsx is NULL.")
  stop_if(!(color_by %in% colnames(scores)), paste0("color_by column not found in metadata: ", color_by))
  p_base <- p_base + aes(color = .data[[color_by]])
}
if (!is.null(shape_by)) {
  stop_if(is.null(metadata_xlsx), "shape_by is set but metadata_xlsx is NULL.")
  stop_if(!(shape_by %in% colnames(scores)), paste0("shape_by column not found in metadata: ", shape_by))
  p_base <- p_base + aes(shape = .data[[shape_by]])
}

ggsave(file.path(outdir, "pca_PC1_PC2.pdf"), p_base, width = 16, height = 12)

# ----------------------------
# Interactive PCA (hover shows sample_id + metadata) WITH same color/shape mappings
# ----------------------------

# Build a tooltip string (works even if color_by/shape_by are NULL)
scores$tooltip <- paste0("sample_id: ", scores$sample_id)

if (!is.null(color_by) && (color_by %in% colnames(scores))) {
  scores$tooltip <- paste0(scores$tooltip, "<br>", color_by, ": ", scores[[color_by]])
}
if (!is.null(shape_by) && (shape_by %in% colnames(scores))) {
  scores$tooltip <- paste0(scores$tooltip, "<br>", shape_by, ": ", scores[[shape_by]])
}

# Start from the SAME base plot (so aesthetics match your PDF)
p_interactive <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_point(aes(text = tooltip), size = 3, alpha = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, title = "PCA: PC1 vs PC2 (interactive)") +
  theme_bw(base_size = 12)

# Add the same color/shape mappings (if set)
if (!is.null(color_by)) {
  stop_if(!(color_by %in% colnames(scores)), paste0("color_by column not found in scores: ", color_by))
  p_interactive <- p_interactive + aes(color = .data[[color_by]])
}
if (!is.null(shape_by)) {
  stop_if(!(shape_by %in% colnames(scores)), paste0("shape_by column not found in scores: ", shape_by))
  p_interactive <- p_interactive + aes(shape = .data[[shape_by]])
}

# Convert to interactive plotly and save HTML
p_int <- ggplotly(p_interactive, tooltip = "text")

htmlwidgets::saveWidget(
  p_int,
  file = file.path(outdir, "pca_PC1_PC2_interactive.html"),
  selfcontained = TRUE
)

message("Wrote: ", file.path(outdir, "pca_PC1_PC2_interactive.html"))

# ----------------------------
# Plot: Loadings barplots (topN for PC1 and PC2)
# ----------------------------
pc1_load <- pca$rotation[, 1]; names(pc1_load) <- rownames(pca$rotation)
pc2_load <- pca$rotation[, 2]; names(pc2_load) <- rownames(pca$rotation)

top1 <- data.frame(feature = names(pc1_load), loading = unname(pc1_load), stringsAsFactors = FALSE)
top1 <- top1[order(abs(top1$loading), decreasing = TRUE), ][seq_len(min(topN, nrow(top1))), ]
top1$PC <- "PC1"

top2 <- data.frame(feature = names(pc2_load), loading = unname(pc2_load), stringsAsFactors = FALSE)
top2 <- top2[order(abs(top2$loading), decreasing = TRUE), ][seq_len(min(topN, nrow(top2))), ]
top2$PC <- "PC2"

top_df <- rbind(top1, top2)

loading_plot <- ggplot(top_df, aes(x = reorder(feature, loading), y = loading)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~PC, scales = "free_y") +
  labs(x = "Feature (id)", y = "Loading",
       title = paste0("Top ", topN, " loadings for PC1 and PC2")) +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "pca_loadings_top.pdf"), loading_plot, width = 10, height = 7)

# ----------------------------
# Plot: PCA biplot with loading arrows (topN features by vector length)
# ----------------------------
rot2 <- pca$rotation[, 1:2, drop = FALSE]
vec_len <- sqrt(rot2[, 1]^2 + rot2[, 2]^2)
ord <- order(vec_len, decreasing = TRUE)
idx <- ord[seq_len(min(topN, length(ord)))]

arrows_df <- data.frame(
  feature = rownames(rot2)[idx],
  PC1 = rot2[idx, 1],
  PC2 = rot2[idx, 2],
  vec_len = vec_len[idx],
  stringsAsFactors = FALSE
)

# Scale arrows to fit score space nicely
score_max1 <- max(abs(scores$PC1), na.rm = TRUE)
score_max2 <- max(abs(scores$PC2), na.rm = TRUE)
load_max1  <- max(abs(arrows_df$PC1), na.rm = TRUE)
load_max2  <- max(abs(arrows_df$PC2), na.rm = TRUE)

scale_factor <- arrow_scale_fraction * min(score_max1 / load_max1, score_max2 / load_max2)

arrows_df$xend <- arrows_df$PC1 * scale_factor
arrows_df$yend <- arrows_df$PC2 * scale_factor

# Label only the strongest vectors (by length in PC1/PC2 loading space)
topN_label <- min(topN_label, nrow(arrows_df))
labels_df <- arrows_df[order(arrows_df$vec_len, decreasing = TRUE), , drop = FALSE]
labels_df <- labels_df[seq_len(topN_label), , drop = FALSE]

# Clean base so arrow layers don't inherit unwanted aesthetics
biplot_base <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.4) +
  geom_point(size = 3, alpha = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, title = paste0("PCA biplot (top ", topN, " loading vectors)")) +
  theme_bw(base_size = 12)

if (!is.null(color_by)) biplot_base <- biplot_base + aes(color = .data[[color_by]])
if (!is.null(shape_by)) biplot_base <- biplot_base + aes(shape = .data[[shape_by]])

biplot_p <- biplot_base +
  geom_segment(
    data = arrows_df,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.4,
    alpha = 0.6
  ) +
  geom_text(
    data = labels_df,
    aes(x = xend, y = yend, label = feature),
    inherit.aes = FALSE,
    size = 3,
    vjust = -0.4,
    check_overlap = TRUE
  )

biplot_path <- file.path(outdir, "pca_biplot_arrows.pdf")
ggsave(biplot_path, biplot_p, width = 16, height = 12)
message("Wrote: ", biplot_path)