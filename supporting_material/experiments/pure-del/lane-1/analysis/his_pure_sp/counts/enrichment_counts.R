# Auto-generated analysis script
suppressPackageStartupMessages({
  library(tidyverse)
  library(GGally)
})

args <- list(
  data_path = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/analysis/his_pure_sp/data.csv",
  samples_path = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/analysis/his_pure_sp/samples.csv",
  save_dir = "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/analysis/his_pure_sp/counts",
  cpm      = FALSE
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

# ---- Load data ----
data <- readr::read_csv(args$data_path, show_col_types = FALSE)
samples <- readr::read_csv(args$samples_path, show_col_types = FALSE)
code_columns <- grep("^code_", colnames(data), value = TRUE)

data = data |>
    dplyr::inner_join(samples, by = "name")

# ---- Optionally compute CPM per library (name) ----
if (isTRUE(args$cpm)) {
  data <- data |>
    dplyr::group_by(name) |>
    dplyr::mutate(count = count / sum(count) * 1e6) |>
    dplyr::ungroup()
}

# ---- Average across replicates ----
data_avg <- data |>
  dplyr::group_by(dplyr::across(dplyr::all_of(code_columns)), group) |>
  dplyr::summarise(mean = mean(count), .groups = "drop")

# ---- Pivot and compute contrasts ----
stats <- data_avg |>
  tidyr::pivot_wider(names_from = group, values_from = mean, values_fill = 0) |>
  dplyr::mutate(
    enrichment = condition - control
    # sticky     = control - naive
  )

# ---- Save outputs ----
readr::write_csv(stats, file.path(args$save_dir, "stats.csv"))

stats |>
  dplyr::arrange(dplyr::desc(enrichment)) |>
  dplyr::slice(1:100) |>
  readr::write_csv(file.path(args$save_dir, "hits.csv"))

# stats |>
#   dplyr::arrange(dplyr::desc(sticky)) |>
#   dplyr::slice(1:100) |>
#   readr::write_csv(file.path(args$save_dir, "sticky.csv"))

# Per-group exports (only those columns present)
present_groups <- intersect(c("condition","control","naive"), colnames(stats))
for (g in present_groups) {
  stats |>
    dplyr::select(dplyr::all_of(code_columns), dplyr::all_of(g)) |>
    dplyr::rename(count = !!rlang::sym(g)) |>
    readr::write_csv(file.path(args$save_dir, paste0(g, ".csv")))
}

for(condition in unique(data$group)){
  g = get_corr_plot(data=data, condition=condition)
  ggsave(file.path(args$save_dir, paste0("correlation_", condition, ".png")), g, width = 8, height = 6)  
}