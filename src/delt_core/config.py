from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

from pydantic import BaseModel, Field, validator


class DemultiplexInitConfig(BaseModel):
    root: Optional[Path] = Field(default=None, description="Root directory to create project structure in")
    experiment_name: str = Field(default="", description="Experiment name")
    selection_file: Path = Field(default=Path("selections/selection.xlsx"))
    library: Path = Field(default=Path("libraries/library.xlsx"))
    fastq_file: Path = Field(default=Path("fastq_files/input.fastq.gz"))
    errors: Optional[List[Tuple[str, float]]] = Field(
        default=None, description="List of (region_name, allowed_error) tuples"
    )


class DemultiplexRunConfig(BaseModel):
    config_file: Path
    write_info_file: bool = False
    write_json_file: bool = True
    fast_dev_run: bool = Field(default=False, description="Run with first 1000 reads only")


class DemultiplexCreateListsConfig(BaseModel):
    config_file: Path
    output_dir: Optional[Path] = None


class DemultiplexComputeCountsConfig(BaseModel):
    config_file: Path
    input_file: Path
    output_dir: Path


class DemultiplexConvertConfig(BaseModel):
    struct_file: Path


class SimulateInitConfig(BaseModel):
    root: Optional[Path] = None
    experiment_name: str = ""
    selection_file: Optional[Path] = None
    library: Optional[Path] = None
    fastq_file: Path = Path("fastq_files/input.fastq.gz")
    output_file: Path = Path("fastq_files/simulation.fastq.gz")
    num_reads: int = 100

    @validator("num_reads")
    def validate_num_reads(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("num_reads must be positive")
        return v


class SimulateRunConfig(BaseModel):
    config_file: Path


class NormalizeRunConfig(BaseModel):
    config_file: Path
    target: str
    control: str


class QcReportConfig(BaseModel):
    experiment_dir: Path


class QcPlotConfig(BaseModel):
    experiment_dir: Path
    output_dir: Optional[Path] = None


class QcCompareLegacyConfig(BaseModel):
    config_file: Path
    legacy_results_dir: Path


class QcAnalyzeCodonsConfig(BaseModel):
    config_file: Path


class ComputeSmilesConfig(BaseModel):
    input_files: List[Path]


class ComputeMergeConfig(BaseModel):
    input_files: List[Path]

    @validator("input_files")
    def require_two_files(cls, v: List[Path]) -> List[Path]:
        if len(v) != 2:
            raise ValueError("merge requires exactly two input files")
        return v


class ComputePropertiesConfig(BaseModel):
    input_file: Path


class ComputePlotConfig(BaseModel):
    input_file: Path


class VizShowConfig(BaseModel):
    path_to_file: Path


class DeltConfig(BaseModel):
    demultiplex_init: Optional[DemultiplexInitConfig] = None
    demultiplex_run: Optional[DemultiplexRunConfig] = None
    demultiplex_create_lists: Optional[DemultiplexCreateListsConfig] = None
    demultiplex_compute_counts: Optional[DemultiplexComputeCountsConfig] = None
    demultiplex_convert: Optional[DemultiplexConvertConfig] = None

    simulate_init: Optional[SimulateInitConfig] = None
    simulate_run: Optional[SimulateRunConfig] = None

    normalize_run: Optional[NormalizeRunConfig] = None

    qc_report: Optional[QcReportConfig] = None
    qc_plot: Optional[QcPlotConfig] = None
    qc_compare_legacy: Optional[QcCompareLegacyConfig] = None
    qc_analyze_codons: Optional[QcAnalyzeCodonsConfig] = None

    compute_smiles: Optional[ComputeSmilesConfig] = None
    compute_merge: Optional[ComputeMergeConfig] = None
    compute_properties: Optional[ComputePropertiesConfig] = None
    compute_plot: Optional[ComputePlotConfig] = None

    viz_show: Optional[VizShowConfig] = None


