from pathlib import Path

import yaml

from delt_hit.demultiplex.preprocess import generate_input_files


def _write_config(tmp_path: Path, input_name: str) -> Path:
    config = {
        "experiment": {
            "save_dir": str(tmp_path),
            "fastq_path": str(tmp_path / input_name),
            "num_cores": 1,
            "name": "demo",
        },
        "structure": [
            {
                "name": "S0",
                "max_error_rate": 0.0,
                "indels": 0,
            }
        ],
        "whitelists": {
            "S0": [{"codon": "ACGT"}],
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path


def test_generate_input_files_preserves_fastq_suffixes_and_stride(tmp_path):
    config_path = _write_config(tmp_path, "sample.fastq.gz")

    script_path = generate_input_files(config_path=config_path)
    script = script_path.read_text()

    assert 'cutadapt_output_files/out.fastq.gz' in script
    assert 'cutadapt_output_files/input.fastq.gz' in script
    assert "sed -n '1~4p'" in script


def test_generate_input_files_preserves_fasta_suffixes_and_stride(tmp_path):
    config_path = _write_config(tmp_path, "sample.fasta.gz")

    script_path = generate_input_files(config_path=config_path)
    script = script_path.read_text()

    assert 'cutadapt_output_files/out.fasta.gz' in script
    assert 'cutadapt_output_files/input.fasta.gz' in script
    assert "sed -n '1~2p'" in script
