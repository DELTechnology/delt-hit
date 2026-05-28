from __future__ import annotations

import math
from pathlib import Path
from textwrap import dedent

import pandas as pd


def compute_zscore_stats(*, counts: pd.DataFrame, library_size: int) -> pd.DataFrame:
    """Compute publication-style z-score statistics from one counts table."""
    if library_size <= 0:
        raise ValueError("Library size must be positive to compute z-score.")

    total_count = float(counts["count"].sum())
    if total_count <= 0:
        raise ValueError("Counts file must have a positive total count.")

    p_exp = 1.0 / float(library_size)
    expected_count = total_count * p_exp
    sigma = math.sqrt(total_count * p_exp * (1.0 - p_exp))
    if sigma == 0:
        raise ValueError("Sigma is zero; cannot compute z-score.")

    stats = counts.copy()
    stats["expected_count"] = expected_count
    stats["sigma"] = sigma
    stats["z_score"] = (stats["count"] - expected_count) / sigma
    stats["norm_z_score"] = stats["z_score"] / math.sqrt(total_count)
    return stats


def zscore_rscript(*, counts_path: Path, library_size: int, save_dir: Path) -> Path:
    """Generate an R script for native z-score analysis."""
    r_script = dedent(
        f"""
        suppressPackageStartupMessages({{
          library(tidyverse)
        }})

        args <- list(
          counts_path = "{counts_path.as_posix()}",
          save_dir = "{save_dir.as_posix()}",
          library_size = {int(library_size)}
        )

        data <- readr::read_tsv(args$counts_path, show_col_types = FALSE)
        code_columns <- grep("^code_", colnames(data), value = TRUE)

        if (!"count" %in% colnames(data)) {{
          stop("Counts file must contain a `count` column.")
        }}
        if (length(code_columns) == 0) {{
          stop("Counts file must contain at least one `code_` column.")
        }}

        selection_total <- sum(data$count)
        if (selection_total <= 0) {{
          stop("Counts file must have a positive total count.")
        }}

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
        """
    ).strip("\n")

    r_path = save_dir / "enrichment_z_score.R"
    r_path.parent.mkdir(parents=True, exist_ok=True)
    r_path.write_text(r_script)
    return r_path
