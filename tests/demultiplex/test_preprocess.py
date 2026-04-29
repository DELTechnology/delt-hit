from pathlib import Path

import yaml

from delt_hit.demultiplex.preprocess import generate_input_files, get_codons


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


def test_get_codons_complements_building_blocks_only():
    whitelists = {
        "B0": [
            {"codon": "GCCTCG", "complement": True},
            {"codon": "AATTCC", "complement": False},
        ],
        "S0": [{"codon": "GCCTCG", "complement": True}],
        "C0": [{"codon": "TTGGCC"}],
    }

    assert get_codons("B0", whitelists) == ["CGGAGC", "AATTCC"]
    assert get_codons("S0", whitelists) == ["GCCTCG"]
    assert get_codons("C0", whitelists) == ["TTGGCC"]


def test_generate_input_files_writes_complemented_building_block_adapters(tmp_path):
    config = {
        "experiment": {
            "save_dir": str(tmp_path),
            "fastq_path": str(tmp_path / "sample.fastq.gz"),
            "num_cores": 1,
            "name": "demo",
        },
        "structure": [
            {
                "name": "B0",
                "max_error_rate": 0.0,
                "indels": 0,
            },
            {
                "name": "S0",
                "max_error_rate": 0.0,
                "indels": 0,
            },
        ],
        "whitelists": {
            "B0": [
                {"index": 0, "codon": "GCCTCG", "complement": True},
                {"index": 1, "codon": "AATTCC", "complement": False},
            ],
            "S0": [{"codon": "GCCTCG"}],
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    generate_input_files(config_path=config_path)

    b0_fastq = (
        tmp_path / "demo" / "demultiplex" / "cutadapt_input_files" / "0-B0.fastq"
    ).read_text()
    s0_fastq = (
        tmp_path / "demo" / "demultiplex" / "cutadapt_input_files" / "1-S0.fastq"
    ).read_text()

    assert b0_fastq == ">0-B0.0\nCGGAGC\n>0-B0.1\nAATTCC"
    assert s0_fastq == ">1-S0.0\nGCCTCG"
