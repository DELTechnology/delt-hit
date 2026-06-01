# Auto-generated analysis script
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(GGally)
})

args <- list(
  data_path    = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-2/analysis/dyna_up/data.csv",
  samples_path = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-2/analysis/dyna_up/samples.csv",
  save_dir     = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-2/analysis/dyna_up/edgeR",
  log          = FALSE
)

# ---- Helper ----
get_corr_plot <- function(data, condition) {
  code_columns <- grep("^code_", colnames(data), value = TRUE)
  pdat <- data %>%
    dplyr::filter(group == condition) %>%
    dplyr::select(dplyr::all_of(code_columns), name, count) %>%
    tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)

  g <- pdat %>%
    dplyr::select(-dplyr::all_of(code_columns)) %>%
    GGally::ggpairs(
      upper = list(continuous = GGally::wrap("smooth", alpha = 0.3, size = 0.2)),
      lower = list(continuous = GGally::wrap("cor", size = 3))
    ) +
    ggplot2::ggtitle(paste("Replicate comparisons for", condition))

  return(g)
}

get_hits.edgeR <- function(data, log = FALSE) {
  code_columns <- grep("^code_", colnames(data), value = TRUE)

  data.wide <- data %>%
    select(-group) %>%
    pivot_wider(names_from = name, values_from = count, values_fill = 0)

  data.row <- data.wide %>%
    select(all_of(code_columns))

  data.col <- data.frame(
    name = data.wide %>% select(-all_of(code_columns), -any_of("id")) %>% colnames(),
    stringsAsFactors = FALSE
  )
  groups <- factor(sapply(data.col$name, get_group_from_name))
  groups <- relevel(groups, "control")
  data.col$group <- groups

  data.counts <- as.matrix(data.wide %>% select(-all_of(code_columns), -any_of("id")))
  rownames(data.counts) <- seq.int(nrow(data.wide))

  y <- DGEList(
    counts  = data.counts,
    samples = data.frame(
      name  = data.col$name,
      group = data.col$group,
      row.names = data.col$name
    )
  )

  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ 0 + y$samples$group)
  colnames(design) <- levels(y$samples$group)

  y   <- estimateDisp(y, design)
  fit <- glmFit(y, design)

  cm <- makeContrasts(
    enrichment = condition - control,
    # sticky     = control - naive,
    levels = design
  )
  lrt.enrichment <- glmLRT(fit, contrast = cm[, 1])
  # lrt.sticky     <- glmLRT(fit, contrast = cm[, 2])

  stats.enrichment <- bind_cols(data.row, lrt.enrichment$table) %>%
    mutate(FDR = p.adjust(PValue, method = "BH"))
  # stats.sticky <- bind_cols(data.row, lrt.sticky$table) %>%
  #   mutate(FDR = p.adjust(PValue, method = "BH"))
  stats.sticky = NULL

  hits.enrichment <- stats.enrichment %>%
    filter(FDR < 0.05, logFC > 0) %>%
    arrange(desc(logFC), FDR)
  # hits.sticky <- stats.sticky %>%
  #   filter(FDR < 0.05, logFC > 0) %>%
  #   arrange(desc(logFC), FDR)
  hits.sticky = NULL

  # Export CPMs (optionally log-transformed)
  counts <- bind_cols(data.row, cpm(y, normalized.lib.sizes = TRUE, log = log, prior.count = 0.5))

  list(
    stats  = list(enrichment = stats.enrichment,
                  sticky     = stats.sticky),
    hits   = list(enrichment = hits.enrichment,
                  sticky     = hits.sticky),
    counts = counts
  )
}

# ---- Load data ----
data    <- read_csv(args$data_path,    show_col_types = FALSE)
samples <- read_csv(args$samples_path, show_col_types = FALSE)
code_columns <- grep("^code_", colnames(data), value = TRUE)

grp_by_name <- setNames(samples$group, samples$name)
get_group_from_name <- function(name) grp_by_name[[name]]

data <- data %>%
  inner_join(samples, by = "name")

# ---- ANALYSIS ----
result.edgeR <- get_hits.edgeR(data = data, log = args$log)

dir.create(args$save_dir, showWarnings = FALSE, recursive = TRUE)

# save stats
for (i in seq_along(result.edgeR$stats)) {
  name  <- names(result.edgeR$stats)[i]
  stats <- result.edgeR$stats[[i]]
  if(is.null(stats)) next
  save.path <- file.path(args$save_dir, paste0(name, "_stats.csv"))
  write_csv(stats, file = save.path)
}

# save hits
for (i in seq_along(result.edgeR$hits)) {
  name <- names(result.edgeR$hits)[i]
  hits <- result.edgeR$hits[[i]]
  if(is.null(hits)) next
  save.path <- file.path(args$save_dir, paste0(name, "_hits.csv"))
  write_csv(hits, file = save.path)
}

# derive selections from result table columns and export per-selection counts
selections <- setdiff(colnames(result.edgeR$counts), code_columns)

# save normalized counts (CPMs; log scale depends on args$log)
for (selection in selections) {
  fname <- paste0(selection, ".csv")
  save.path <- file.path(args$save_dir, fname)
  result.edgeR$counts %>%
    select(all_of(code_columns), all_of(selection)) %>%
    mutate(count = .data[[selection]]) %>%
    select(-all_of(selection)) %>%
    write_csv(file = save.path)
}

counts.norm <- result.edgeR$counts %>%
  pivot_longer(
    cols = -all_of(code_columns),
    names_to = "name",
    values_to = "count"
  )
counts.norm$group = sapply(counts.norm$name, get_group_from_name)
for(condition in unique(counts.norm$group)){
  g = get_corr_plot(data=counts.norm, condition=condition)
  ggsave(file.path(args$save_dir, paste0("correlation_", condition, ".png")), g, width = 8, height = 6)  
}