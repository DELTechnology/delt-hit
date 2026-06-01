from __future__ import annotations

from pathlib import Path
from textwrap import dedent

from loguru import logger


def edgeR_rscript(*, data_path: Path, samples_path: Path, log: bool = False, save_dir: Path):
    """Generate an edgeR analysis R script."""
    log_flag = "TRUE" if log else "FALSE"

    r_script = dedent(f"""
        # Auto-generated analysis script
        suppressPackageStartupMessages({{
          library(tidyverse)
          library(edgeR)
          library(limma)
          library(GGally)
        }})

        args <- list(
          data_path    = "{data_path.as_posix()}",
          samples_path = "{samples_path.as_posix()}",
          save_dir     = "{save_dir.as_posix()}",
          log          = {log_flag}
        )

        # ---- Helper ----
        get_corr_plot <- function(data, condition) {{
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
        }}

        get_hits.edgeR <- function(data, log = FALSE) {{
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
            levels = design
          )
          lrt.enrichment <- glmLRT(fit, contrast = cm[, 1])

          stats.enrichment <- bind_cols(data.row, lrt.enrichment$table) %>%
            mutate(FDR = p.adjust(PValue, method = "BH"))
          hits.enrichment <- stats.enrichment %>%
            filter(FDR < 0.05, logFC > 0) %>%
            arrange(desc(logFC), FDR)

          counts <- bind_cols(data.row, cpm(y, normalized.lib.sizes = TRUE, log = log, prior.count = 0.5))

          list(
            stats  = list(enrichment = stats.enrichment),
            hits   = list(enrichment = hits.enrichment),
            counts = counts
          )
        }}

        data    <- read_csv(args$data_path,    show_col_types = FALSE)
        samples <- read_csv(args$samples_path, show_col_types = FALSE)
        code_columns <- grep("^code_", colnames(data), value = TRUE)

        grp_by_name <- setNames(samples$group, samples$name)
        get_group_from_name <- function(name) grp_by_name[[name]]

        data <- data %>%
          inner_join(samples, by = "name")

        result.edgeR <- get_hits.edgeR(data = data, log = args$log)

        dir.create(args$save_dir, showWarnings = FALSE, recursive = TRUE)

        stats <- result.edgeR$stats[["enrichment"]]
        if (!is.null(stats)) {{
          write_csv(stats, file.path(args$save_dir, "stats.csv"))
        }}

        hits <- result.edgeR$hits[["enrichment"]]
        if (!is.null(hits)) {{
          write_csv(hits, file.path(args$save_dir, "hits.csv"))
        }}

        selections <- setdiff(colnames(result.edgeR$counts), code_columns)

        for (selection in selections) {{
          fname <- paste0(selection, ".csv")
          save.path <- file.path(args$save_dir, fname)
          result.edgeR$counts %>%
            select(all_of(code_columns), all_of(selection)) %>%
            mutate(count = .data[[selection]]) %>%
            select(-all_of(selection)) %>%
            write_csv(file = save.path)
        }}

        counts.norm <- result.edgeR$counts %>%
          pivot_longer(
            cols = -all_of(code_columns),
            names_to = "name",
            values_to = "count"
          )
        counts.norm$group = sapply(counts.norm$name, get_group_from_name)
        for(condition in unique(counts.norm$group)){{
          g = get_corr_plot(data=counts.norm, condition=condition)
          ggsave(file.path(args$save_dir, paste0("correlation_", condition, ".png")), g, width = 8, height = 6)
        }}
    """).strip("\n")

    r_path = save_dir / "enrichment_edgeR.R"
    r_path.parent.mkdir(parents=True, exist_ok=True)
    r_path.write_text(r_script)
    logger.info(f"Created R script at {r_path}")


def counts_rscript(*, data_path: Path, samples_path: Path, cpm, save_dir: Path):
    """Generate a simple counts-based analysis R script."""
    r_cpm_flag = "TRUE" if cpm else "FALSE"

    r_script = dedent(f"""
        # Auto-generated analysis script
        suppressPackageStartupMessages({{
          library(tidyverse)
          library(GGally)
        }})

        args <- list(
          data_path = "{data_path.as_posix()}",
          samples_path = "{samples_path.as_posix()}",
          save_dir = "{save_dir.as_posix()}",
          cpm      = {r_cpm_flag}
        )

        get_corr_plot <- function(data, condition) {{
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
        }}

        data <- readr::read_csv(args$data_path, show_col_types = FALSE)
        samples <- readr::read_csv(args$samples_path, show_col_types = FALSE)
        code_columns <- grep("^code_", colnames(data), value = TRUE)

        data = data |>
            dplyr::inner_join(samples, by = "name")

        if (isTRUE(args$cpm)) {{
          data <- data |>
            dplyr::group_by(name) |>
            dplyr::mutate(count = count / sum(count) * 1e6) |>
            dplyr::ungroup()
        }}

        data_avg <- data |>
          dplyr::group_by(dplyr::across(dplyr::all_of(code_columns)), group) |>
          dplyr::summarise(mean = mean(count), .groups = "drop")

        stats <- data_avg |>
          tidyr::pivot_wider(names_from = group, values_from = mean, values_fill = 0) |>
          dplyr::mutate(
            enrichment = condition - control
          )

        readr::write_csv(stats, file.path(args$save_dir, "stats.csv"))

        stats |>
          dplyr::arrange(dplyr::desc(enrichment)) |>
          dplyr::slice(1:100) |>
          readr::write_csv(file.path(args$save_dir, "hits.csv"))

        present_groups <- intersect(c("condition","control","naive"), colnames(stats))
        for (g in present_groups) {{
          stats |>
            dplyr::select(dplyr::all_of(code_columns), dplyr::all_of(g)) |>
            dplyr::rename(count = !!rlang::sym(g)) |>
            readr::write_csv(file.path(args$save_dir, paste0(g, ".csv")))
        }}

        for(condition in unique(data$group)){{
          g = get_corr_plot(data=data, condition=condition)
          ggsave(file.path(args$save_dir, paste0("correlation_", condition, ".png")), g, width = 8, height = 6)
        }}
    """).strip("\n")

    r_path = save_dir / "enrichment_counts.R"
    r_path.parent.mkdir(parents=True, exist_ok=True)
    r_path.write_text(r_script)
    logger.info(f"Created R script at {r_path}")
