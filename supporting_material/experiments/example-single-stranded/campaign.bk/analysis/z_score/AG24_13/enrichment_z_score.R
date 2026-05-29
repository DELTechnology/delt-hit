suppressPackageStartupMessages({
  library(tidyverse)
})

args <- list(
  counts_path = "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/example-single-stranded/campaign/selections/AG24_13/counts.txt",
  save_dir = "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/example-single-stranded/campaign/analysis/z_score/AG24_13",
  library_size = 168337
)

data <- readr::read_tsv(args$counts_path, show_col_types = FALSE)
code_columns <- grep("^code_", colnames(data), value = TRUE)

if (!"count" %in% colnames(data)) {
  stop("Counts file must contain a `count` column.")
}
if (length(code_columns) == 0) {
  stop("Counts file must contain at least one `code_` column.")
}

selection_total <- sum(data$count)
if (selection_total <= 0) {
  stop("Counts file must have a positive total count.")
}

p_exp <- 1 / args$library_size
expected_count <- selection_total * p_exp
sigma <- sqrt(selection_total * p_exp * (1 - p_exp))

stats <- data %>%
  mutate(
    expected_count = expected_count,
    sigma = sigma,
    z_score = (count - expected_count) / sigma,
    norm_z_score = z_score / sqrt(selection_total)
  )

dir.create(args$save_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(stats, file.path(args$save_dir, "stats.csv"))
stats %>%
  arrange(desc(norm_z_score)) %>%
  slice(1:100) %>%
  readr::write_csv(file.path(args$save_dir, "hits.csv"))