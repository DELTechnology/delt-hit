from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


def _write_counts(path: Path, rows: list[tuple[str, str, int]]) -> None:
    frame = pd.DataFrame(rows, columns=["code_1", "code_2", "count"])
    frame["id"] = [f"cmpd_{idx}" for idx in range(len(frame))]
    frame.to_csv(path, sep="\t", index=False)


def test_delt_hit_analyse_enrichment_zscore_cli(tmp_path: Path):
    protein_counts = tmp_path / "protein.tsv"
    control_counts = tmp_path / "control.tsv"
    _write_counts(protein_counts, [("A", "A", 90), ("A", "B", 10)])
    _write_counts(control_counts, [("A", "A", 10), ("A", "B", 90)])

    output_root = tmp_path / "out"
    config = tmp_path / "analysis.yaml"
    config.write_text(
        f"""
experiments:
  - name: demo
    save_dir: {output_root}
    library_size: 4
    selections:
      - name: protein_1
        counts_path: {protein_counts}
        group: protein
      - name: control_1
        counts_path: {control_counts}
        group: no_protein
""".strip()
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "delt_hit.cli.main",
            "analyse",
            "enrichment",
            "--config_path",
            str(config),
            "--name",
            "demo",
            "--method",
            "zscore",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (output_root / "demo" / "zscore" / "stats.csv").exists()
    assert (output_root / "demo" / "zscore" / "hits.csv").exists()
