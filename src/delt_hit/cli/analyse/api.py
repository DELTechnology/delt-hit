from pathlib import Path

from loguru import logger

from delt_hit.analyse.config import derive_library_size, load_analysis_experiment, load_project_config
from delt_hit.analyse.data_prep import load_zscore_counts, prepare_replicate_analysis_data
from delt_hit.analyse.scripts import counts_rscript, edgeR_rscript
from delt_hit.analyse.zscore import zscore_rscript

class Analyse:
    def enrichment(
        self,
        *,
        method: str = "counts",
        save_dir: Path | None = None,
        analysis_config: Path | None = None,
        name: str | None = None,
        config_path: Path | None = None,
        counts: Path | None = None,
    ):
        """Generate enrichment analysis scripts for an experiment.

        Args:
            method: Analysis method.
            save_dir: Explicit output directory for method artifacts.
            analysis_config: Analysis YAML for replicate-based methods.
            name: Experiment name for replicate-based methods.
            config_path: Init-generated project config for z-score.
            counts: Counts table path for z-score.
        """
        self._validate_enrichment_inputs(
            method=method,
            save_dir=save_dir,
            analysis_config=analysis_config,
            name=name,
            config_path=config_path,
            counts=counts,
        )
        assert save_dir is not None
        save_dir = Path(save_dir).expanduser().resolve()

        if method in {"counts", "edgeR", "DESeq2"}:
            assert analysis_config is not None
            assert name is not None
            method_dir = save_dir / method / name
            exp = load_analysis_experiment(analysis_config=analysis_config, name=name)
            data_path = method_dir / "data.csv"
            samples_path = method_dir / "samples.csv"
            prepare_replicate_analysis_data(exp=exp, data_path=data_path, samples_path=samples_path)
            logger.info(f"Prepared data at {data_path} and samples at {samples_path}")

            match method:
                case "counts":
                    counts_rscript(data_path=data_path, samples_path=samples_path, cpm=False, save_dir=method_dir)
                case "edgeR":
                    edgeR_rscript(data_path=data_path, samples_path=samples_path, log=False, save_dir=method_dir)
                case "DESeq2":
                    raise NotImplementedError("DESeq2 analysis is not implemented.")
            return

        assert config_path is not None
        assert counts is not None
        assert name is not None
        project_config = load_project_config(config_path=config_path)
        library_size = derive_library_size(project_config=project_config)
        load_zscore_counts(counts_path=counts)

        method_dir = save_dir / "z_score" / name
        method_dir.mkdir(parents=True, exist_ok=True)
        script_path = zscore_rscript(
            counts_path=Path(counts).expanduser().resolve(),
            library_size=library_size,
            save_dir=method_dir,
        )
        logger.info(f"Created z-score analysis script at {script_path}")

    @staticmethod
    def _validate_enrichment_inputs(
        *,
        method: str,
        save_dir: Path | None,
        analysis_config: Path | None,
        name: str | None,
        config_path: Path | None,
        counts: Path | None,
    ) -> None:
        if save_dir is None:
            raise ValueError("`save_dir` is required for all analysis methods.")

        if method in {"counts", "edgeR", "DESeq2"}:
            if analysis_config is None:
                raise ValueError(f"`analysis_config` is required for method `{method}`.")
            if name is None:
                raise ValueError(f"`name` is required for method `{method}`.")
            if config_path is not None:
                raise ValueError(f"`config_path` is not supported for method `{method}`.")
            if counts is not None:
                raise ValueError(f"`counts` is not supported for method `{method}`.")
            return

        if method == "z_score":
            if config_path is None:
                raise ValueError("`config_path` is required for method `z_score`.")
            if counts is None:
                raise ValueError("`counts` is required for method `z_score`.")
            if name is None:
                raise ValueError("`name` is required for method `z_score`.")
            if analysis_config is not None:
                raise ValueError("`analysis_config` is not supported for method `z_score`.")
            return

        raise ValueError(f"Unsupported analysis method: {method}")
