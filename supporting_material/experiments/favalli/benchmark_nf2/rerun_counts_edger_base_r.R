#!/usr/bin/env Rscript
# Recompute DELT-Hit simple-count and edgeR outputs without tidyverse/GGally.
#
# Usage from this benchmark_nf2 directory after an Analyse.enrichment run has
# created outputs/<experiment>/data.csv and outputs/<experiment>/samples.csv:
#
#   Rscript rerun_counts_edger_base_r.R outputs/ss_caix
#
# The generated DELT-Hit R scripts use tidyverse/GGally for data reshaping and
# optional correlation plots. This compatibility script keeps the statistical
# edgeR model unchanged while using base R for I/O and reshaping.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "outputs/ss_caix"
base_dir <- normalizePath(base_dir, mustWork = TRUE)

data_path <- file.path(base_dir, "data.csv")
samples_path <- file.path(base_dir, "samples.csv")
counts_dir <- file.path(base_dir, "counts")
edger_dir <- file.path(base_dir, "edgeR")

if (!file.exists(data_path)) {
  stop("Missing data.csv: ", data_path)
}
if (!file.exists(samples_path)) {
  stop("Missing samples.csv: ", samples_path)
}

dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(edger_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.csv(samples_path, stringsAsFactors = FALSE, check.names = FALSE)
code_columns <- grep("^code_", names(data), value = TRUE)
if (length(code_columns) == 0) {
  stop("No code_* columns found in ", data_path)
}

data <- merge(data, samples, by = "name", all.x = FALSE, all.y = FALSE, sort = FALSE)

# Simple-count workflow: mean raw count per group, then protein - no_protein.
data_avg <- aggregate(data$count, by = data[c(code_columns, "group")], FUN = mean)
names(data_avg)[ncol(data_avg)] <- "mean"

stats_counts <- reshape(
  data_avg,
  idvar = code_columns,
  timevar = "group",
  direction = "wide"
)
names(stats_counts) <- sub("^mean\\.", "", names(stats_counts))

for (group_name in c("protein", "no_protein", "naive")) {
  if (!group_name %in% names(stats_counts)) {
    stats_counts[[group_name]] <- 0
  }
}
stats_counts[is.na(stats_counts)] <- 0
stats_counts$enrichment <- stats_counts$protein - stats_counts$no_protein
stats_counts <- stats_counts[order(-stats_counts$enrichment), ]

write.csv(stats_counts, file.path(counts_dir, "stats.csv"), row.names = FALSE)
write.csv(head(stats_counts, 100), file.path(counts_dir, "hits.csv"), row.names = FALSE)

for (group_name in intersect(c("protein", "no_protein", "naive"), names(stats_counts))) {
  group_table <- stats_counts[c(code_columns, group_name)]
  names(group_table)[ncol(group_table)] <- "count"
  write.csv(group_table, file.path(counts_dir, paste0(group_name, ".csv")), row.names = FALSE)
}

# edgeR workflow: wide sparse count matrix with absent counts set to zero.
wide <- reshape(
  data[c(code_columns, "id", "name", "count")],
  idvar = c(code_columns, "id"),
  timevar = "name",
  direction = "wide"
)
selection_cols <- grep("^count\\.", names(wide), value = TRUE)
colnames(wide)[match(selection_cols, names(wide))] <- sub("^count\\.", "", selection_cols)

selection_names <- setdiff(names(wide), c(code_columns, "id"))
wide[selection_names][is.na(wide[selection_names])] <- 0

sample_groups <- samples$group[match(selection_names, samples$name)]
sample_groups <- factor(sample_groups, levels = c("no_protein", "protein"))

counts_mat <- as.matrix(wide[selection_names])
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- seq_len(nrow(counts_mat))

y <- DGEList(
  counts = counts_mat,
  samples = data.frame(
    name = selection_names,
    group = sample_groups,
    row.names = selection_names
  )
)
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + y$samples$group)
colnames(design) <- levels(y$samples$group)

y <- estimateDisp(y, design)
fit <- glmFit(y, design)
contrast_matrix <- makeContrasts(enrichment = protein - no_protein, levels = design)
lrt <- glmLRT(fit, contrast = contrast_matrix[, 1])

stats_edger <- cbind(wide[code_columns], lrt$table)
stats_edger$FDR <- p.adjust(stats_edger$PValue, method = "BH")
stats_edger <- stats_edger[order(-stats_edger$logFC, stats_edger$FDR), ]

write.csv(stats_edger, file.path(edger_dir, "enrichment_stats.csv"), row.names = FALSE)

hits_edger <- stats_edger[stats_edger$FDR < 0.05 & stats_edger$logFC > 0, ]
write.csv(hits_edger, file.path(edger_dir, "enrichment_hits.csv"), row.names = FALSE)

cpm_counts <- cbind(
  wide[code_columns],
  cpm(y, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0.5)
)
for (selection_name in selection_names) {
  selection_table <- cpm_counts[c(code_columns, selection_name)]
  names(selection_table)[ncol(selection_table)] <- "count"
  write.csv(selection_table, file.path(edger_dir, paste0(selection_name, ".csv")), row.names = FALSE)
}

cat("Wrote counts outputs to ", counts_dir, "\n", sep = "")
cat("Wrote edgeR outputs to ", edger_dir, "\n", sep = "")
cat("edgeR common dispersion: ", y$common.dispersion, "\n", sep = "")
