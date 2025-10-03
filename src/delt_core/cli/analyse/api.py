from pathlib import Path
from textwrap import dedent

import pandas as pd

from delt_core.utils import read_yaml, write_yaml


# config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-3/config.yaml')
# cfg = read_yaml(config_path)

class Analyse:

    def prepare(self, config_path: Path, name: str):
        cfg = read_yaml(config_path)
        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        selections_dir = exp_dir / 'selections'

        analysis_dir = exp_dir / 'analyses' / name

        data_path = analysis_dir / 'data.csv'
        samples_path = analysis_dir / 'samples.csv'

        prepare_data(cfg=cfg, name=name, selections_dir=selections_dir,
                     data_path=data_path, samples_path=samples_path)

        return data_path, samples_path, analysis_dir

    def enrichment(self, *, config_path: Path, name: str, method: str = 'counts'):
        data_path, samples_path, analysis_dir = self.prepare(config_path=config_path, name=name)

        match method:
            case 'counts':
                counts_rscript(data_path=data_path,
                               samples_path=samples_path,
                               cpm=False,
                               save_dir=analysis_dir / 'counts')
            case 'edgeR':
                edgeR_rscript(data_path=data_path,
                              samples_path=samples_path,
                              log=False,
                              save_dir=analysis_dir / 'edgeR')
            case 'DESeq2':
                pass

    def run(self, config_path: Path):
        pass


def correlation_rscript(*, data_path: Path, samples_path: Path, cpm, save_dir: Path):
    pass


def edgeR_rscript(*, data_path: Path, samples_path: Path, log: bool = False, save_dir: Path):
    from textwrap import dedent

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
          pdat <- data %>%
            dplyr::filter(group == condition) %>%
            dplyr::select(code_1, code_2, name, count) %>%
            tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)
        
          g <- pdat %>%
            dplyr::select(-code_1, -code_2) %>%
            GGally::ggpairs(
              upper = list(continuous = GGally::wrap("smooth", alpha = 0.3, size = 0.2)),
              lower = list(continuous = GGally::wrap("cor", size = 3))
            ) +
            ggplot2::ggtitle(paste("Replicate comparisons for", condition))
        
          return(g)
        }}
        
        get_hits.edgeR <- function(data, log = FALSE) {{

          data.wide <- data %>%
            select(-group) %>%
            pivot_wider(names_from = name, values_from = count, values_fill = 0)

          data.row <- data.wide %>%
            select(code_1, code_2)

          data.col <- data.frame(
            name = data.wide %>% select(-code_1, -code_2) %>% colnames(),
            stringsAsFactors = FALSE
          )
          groups <- factor(sapply(data.col$name, get_group_from_name))
          groups <- relevel(groups, "no_protein")
          data.col$group <- groups

          data.counts <- as.matrix(data.wide %>% select(-code_1, -code_2))
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
            enrichment = protein - no_protein,
            # sticky     = no_protein - naive,
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
        }}

        # ---- Load data ----
        data    <- read_csv(args$data_path,    show_col_types = FALSE)
        samples <- read_csv(args$samples_path, show_col_types = FALSE)

        grp_by_name <- setNames(samples$group, samples$name)
        get_group_from_name <- function(name) grp_by_name[[name]]

        data <- data %>%
          inner_join(samples, by = "name")

        # ---- ANALYSIS ----
        result.edgeR <- get_hits.edgeR(data = data, log = args$log)

        dir.create(args$save_dir, showWarnings = FALSE, recursive = TRUE)

        # save stats
        for (i in seq_along(result.edgeR$stats)) {{
          name  <- names(result.edgeR$stats)[i]
          stats <- result.edgeR$stats[[i]]
          if(is.null(stats)) next
          save.path <- file.path(args$save_dir, paste0(name, "_stats.csv"))
          write_csv(stats, file = save.path)
        }}

        # save hits
        for (i in seq_along(result.edgeR$hits)) {{
          name <- names(result.edgeR$hits)[i]
          hits <- result.edgeR$hits[[i]]
          if(is.null(hits)) next
          save.path <- file.path(args$save_dir, paste0(name, "_hits.csv"))
          write_csv(hits, file = save.path)
        }}

        # derive selections from result table columns and export per-selection counts
        selections <- setdiff(colnames(result.edgeR$counts), c("code_1", "code_2"))

        # save normalized counts (CPMs; log scale depends on args$log)
        for (selection in selections) {{
          fname <- paste0(selection, ".csv")
          save.path <- file.path(args$save_dir, fname)
          result.edgeR$counts %>%
            select(code_1, code_2, all_of(selection)) %>%
            mutate(count = .data[[selection]]) %>%
            select(-all_of(selection)) %>%
            write_csv(file = save.path)
        }}
        
        counts.norm <- result.edgeR$counts %>%
          pivot_longer(
            cols = -c(code_1, code_2),
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


def counts_rscript(*, data_path: Path, samples_path: Path, cpm, save_dir: Path):
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
        
        # ---- Helper ----
        get_corr_plot <- function(data, condition) {{
          pdat <- data %>%
            dplyr::filter(group == condition) %>%
            dplyr::select(code_1, code_2, name, count) %>%
            tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)
        
          g <- pdat %>%
            dplyr::select(-code_1, -code_2) %>%
            GGally::ggpairs(
              upper = list(continuous = GGally::wrap("smooth", alpha = 0.3, size = 0.2)),
              lower = list(continuous = GGally::wrap("cor", size = 3))
            ) +
            ggplot2::ggtitle(paste("Replicate comparisons for", condition))
        
          return(g)
        }}

        # ---- Load data ----
        data <- readr::read_csv(args$data_path, show_col_types = FALSE)
        samples <- readr::read_csv(args$samples_path, show_col_types = FALSE)
        
        data = data |>
            dplyr::inner_join(samples, by = "name")

        # ---- Optionally compute CPM per library (name) ----
        if (isTRUE(args$cpm)) {{
          data <- data |>
            dplyr::group_by(name) |>
            dplyr::mutate(count = count / sum(count) * 1e6) |>
            dplyr::ungroup()
        }}

        # ---- Average across replicates ----
        data_avg <- data |>
          dplyr::group_by(code_1, code_2, group) |>
          dplyr::summarise(mean = mean(count), .groups = "drop")

        # ---- Pivot and compute contrasts ----
        stats <- data_avg |>
          tidyr::pivot_wider(names_from = group, values_from = mean, values_fill = 0) |>
          dplyr::mutate(
            enrichment = protein - no_protein
            # sticky     = no_protein - naive
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
        present_groups <- intersect(c("protein","no_protein","naive"), colnames(stats))
        for (g in present_groups) {{
          stats |>
            dplyr::select(code_1, code_2, dplyr::all_of(g)) |>
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


def prepare_data(cfg, name: str, selections_dir: Path, data_path: Path, samples_path: Path):
    selections = cfg['analyses'][name]['selections']
    samples = []
    data = []
    for sel in selections:
        meta = cfg['selections'][sel]
        meta['name'] = sel
        samples.append(meta)
        counts = pd.read_csv(selections_dir / sel / "counts.txt", delimiter='\t')
        counts['name'] = sel
        data.append(counts)

    samples = pd.DataFrame(samples)[['name', 'group']]
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(samples_path, index=False)

    data = pd.concat(data, axis=0)
    data_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(data_path, index=False)
