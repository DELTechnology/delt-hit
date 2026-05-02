import json
from pathlib import Path

from delt_hit.quality_control.plot_codon_hits import plot_hits


def _write_minimal_cutadapt_report(path: Path) -> None:
    report = {
        "adapters_read1": [
            {
                "name": "S0.0",
                "five_prime_end": {
                    "trimmed_lengths": [
                        {"counts": {"0": 10, "1": 2}},
                    ]
                },
            }
        ],
        "read_counts": {
            "input": 20,
            "output": 12,
        },
    }
    path.write_text(json.dumps(report))


def test_plot_hits_writes_pdf_png_and_parquet(tmp_path: Path) -> None:
    output_dir = tmp_path / "cutadapt_output_files"
    save_dir = tmp_path / "qc"
    output_dir.mkdir()
    save_dir.mkdir()

    _write_minimal_cutadapt_report(output_dir / "0-S0.cutadapt.json")

    plot_hits(output_dir=output_dir, save_dir=save_dir)

    assert (save_dir / "hits_S0.parquet").exists()
    assert (save_dir / "hits_S0.pdf").exists()
    assert (save_dir / "hits_S0.png").exists()
