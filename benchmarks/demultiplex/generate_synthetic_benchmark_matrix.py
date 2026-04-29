#!/usr/bin/env python3
"""Generate one or more synthetic FASTQ benchmark datasets from named presets."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BENCHMARKS_ROOT = Path(__file__).resolve().parent
if str(BENCHMARKS_ROOT) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS_ROOT))

from generate_synthetic_fastq import main as generate_dataset


DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "benchmarks" / "demultiplex" / "data"
DEFAULT_PROFILE = "canonical"
DEFAULT_NUM_ERRORS = 0
DEPTH_LABELS = ("1m", "10m", "100m", "1000m")


@dataclass(frozen=True)
class DatasetSpec:
    experiment_name: str
    num_cycles: int
    building_blocks_per_cycle: int
    num_reads_per_compound: int

    @property
    def expected_compounds(self) -> int:
        return self.building_blocks_per_cycle ** self.num_cycles

    @property
    def expected_reads(self) -> int:
        return self.expected_compounds * self.num_reads_per_compound


def make_dataset_spec(*, num_cycles: int, building_blocks_per_cycle: int, depth_label: str) -> DatasetSpec:
    total_reads_by_depth = {
        "1m": 1_000_000,
        "10m": 10_000_000,
        "100m": 100_000_000,
        "1000m": 1_000_000_000,
    }
    try:
        total_reads = total_reads_by_depth[depth_label]
    except KeyError as exc:
        raise ValueError(f"Unsupported depth label: {depth_label}") from exc

    expected_compounds = building_blocks_per_cycle ** num_cycles
    if total_reads % expected_compounds != 0:
        raise ValueError(
            f"Depth {depth_label} is not divisible by the number of compounds "
            f"({expected_compounds}) for {num_cycles} cycles x {building_blocks_per_cycle} BBs"
        )

    return DatasetSpec(
        experiment_name=f"synthetic_{num_cycles}cycle_{building_blocks_per_cycle}bbpc_{depth_label}",
        num_cycles=num_cycles,
        building_blocks_per_cycle=building_blocks_per_cycle,
        num_reads_per_compound=total_reads // expected_compounds,
    )


def canonical_specs() -> list[DatasetSpec]:
    return [
        make_dataset_spec(
            num_cycles=num_cycles,
            building_blocks_per_cycle=10,
            depth_label=depth_label,
        )
        for num_cycles in (2, 3, 4)
        for depth_label in DEPTH_LABELS
    ]


def two_cycle_large_library_specs() -> list[DatasetSpec]:
    return [
        make_dataset_spec(
            num_cycles=2,
            building_blocks_per_cycle=building_blocks_per_cycle,
            depth_label=depth_label,
        )
        for building_blocks_per_cycle in (100, 1000)
        for depth_label in DEPTH_LABELS
    ]


PROFILES: dict[str, list[DatasetSpec]] = {
    "canonical": canonical_specs(),
    "two_cycle_large_libraries": two_cycle_large_library_specs(),
    "all": canonical_specs() + two_cycle_large_library_specs(),
}


def get_profile_specs(profile: str) -> list[DatasetSpec]:
    try:
        return PROFILES[profile]
    except KeyError as exc:
        raise ValueError(f"Unknown profile: {profile}") from exc


def get_profile_dataset_names(profile: str) -> list[str]:
    return [spec.experiment_name for spec in get_profile_specs(profile)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=tuple(PROFILES),
        default=DEFAULT_PROFILE,
        help="Named dataset matrix to generate.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory that will contain one folder per generated dataset.",
    )
    parser.add_argument(
        "--num-errors",
        type=int,
        choices=(0, 1),
        default=DEFAULT_NUM_ERRORS,
        help="Number of substitution errors to inject into each barcode region.",
    )
    parser.add_argument(
        "--errors-in-bb-only",
        action="store_true",
        help="When set with --num-errors 1, mutate only building-block tags.",
    )
    parser.add_argument(
        "--compressed",
        default="true",
        choices=("true", "false"),
        help="Whether to gzip compress FASTQ output.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the dataset plan without writing any files.",
    )
    return parser.parse_args()


def print_plan(specs: list[DatasetSpec]) -> None:
    for spec in specs:
        print(
            f"{spec.experiment_name}: "
            f"{spec.num_cycles} cycles, "
            f"{spec.building_blocks_per_cycle} BB/cycle, "
            f"{spec.expected_compounds} compounds, "
            f"{spec.num_reads_per_compound} reads/compound, "
            f"{spec.expected_reads} total reads"
        )


def main() -> None:
    args = parse_args()
    specs = get_profile_specs(args.profile)
    print_plan(specs)
    if args.dry_run:
        return

    for spec in specs:
        generate_dataset(
            num_cycles=spec.num_cycles,
            building_blocks_per_cycle=spec.building_blocks_per_cycle,
            num_reads_per_compound=spec.num_reads_per_compound,
            num_errors=args.num_errors,
            errors_in_bb_only=args.errors_in_bb_only,
            output_dir=args.output_dir,
            experiment_name=spec.experiment_name,
            compressed=args.compressed == "true",
        )


if __name__ == "__main__":
    main()
