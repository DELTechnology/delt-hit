# Auto-generated analysis script
suppressPackageStartupMessages({
  library(tidyverse)
  library(GGally)
})

args <- list(
  data_path = "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-1-fasta/analysis/counts/ca9_ss/data.csv",
  samples_path = "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-1-fasta/analysis/counts/ca9_ss/samples.csv",
  save_dir = "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-1-fasta/analysis/counts/ca9_ss",
  cpm      = FALSE
)

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

data <- readr::read_csv(args$data_path, show_col_types = FALSE)
samples <- readr::read_csv(args$samples_path, show_col_types = FALSE)
code_columns <- grep("^code_", colnames(data), value = TRUE)

data = data |>
    dplyr::inner_join(samples, by = "name")

if (isTRUE(args$cpm)) {
  data <- data |>
    dplyr::group_by(name) |>
    dplyr::mutate(count = count / sum(count) * 1e6) |>
    dplyr::ungroup()
}

data_avg <- data |>
  dplyr::group_by(dplyr::across(dplyr::all_of(code_columns)), group) |>
  dplyr::summarise(mean = mean(count), .groups = "drop")

stats <- data_avg |>
  tidyr::pivot_wider(names_from = group, values_from = mean, values_fill = 0) |>
  dplyr::mutate(
    enrichment = protein - no_protein
  )

readr::write_csv(stats, file.path(args$save_dir, "stats.csv"))

stats |>
  dplyr::arrange(dplyr::desc(enrichment)) |>
  dplyr::slice(1:100) |>
  readr::write_csv(file.path(args$save_dir, "hits.csv"))

present_groups <- intersect(c("protein","no_protein","naive"), colnames(stats))
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