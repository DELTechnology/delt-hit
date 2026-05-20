from pathlib import Path

import yaml

from delt_hit.demultiplex.preprocess import generate_input_files, get_codons
from delt_hit.demultiplex.validation import Region


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
    assert 'gzip -cd ' in script
    assert "awk 'NR % 4 == 1'" in script
    assert "sed -n '1~4p'" not in script


def test_generate_input_files_preserves_fasta_suffixes_and_stride(tmp_path):
    config_path = _write_config(tmp_path, "sample.fasta.gz")

    script_path = generate_input_files(config_path=config_path)
    script = script_path.read_text()

    assert 'cutadapt_output_files/out.fasta.gz' in script
    assert 'cutadapt_output_files/input.fasta.gz' in script
    assert 'gzip -cd ' in script
    assert "awk 'NR % 2 == 1'" in script
    assert "sed -n '1~2p'" not in script


def test_generate_input_files_normalizes_dotted_fastq_name(tmp_path):
    config_path = _write_config(tmp_path, "20190812.A-1907_NF2GB2_s1_R1.fastq.gz")

    script_path = generate_input_files(config_path=config_path)
    script = script_path.read_text()

    assert 'cutadapt_output_files/out.fastq.gz' in script
    assert 'cutadapt_output_files/input.fastq.gz' in script
    assert 'out.A-1907_NF2GB2_s1_R1.fastq.gz' not in script
    assert 'input.A-1907_NF2GB2_s1_R1.fastq.gz' not in script


def test_generate_input_files_normalizes_short_fasta_extension(tmp_path):
    config_path = _write_config(tmp_path, "sample.fa.gz")

    script_path = generate_input_files(config_path=config_path)
    script = script_path.read_text()

    assert 'cutadapt_output_files/out.fasta.gz' in script
    assert 'cutadapt_output_files/input.fasta.gz' in script
    assert "awk 'NR % 2 == 1'" in script


def test_get_codons_applies_region_level_complement_to_all_codons():
    region = Region(
        name="B0",
        index=0,
        codons=[],
        max_error_rate=0.0,
        indels=0,
        complement=True,
    )
    whitelists = {
        "B0": [
            {"codon": "GCCTCG"},
            {"codon": "AATTCC"},
        ],
        "C0": [{"codon": "TTGGCC"}],
    }

    assert get_codons(region, whitelists) == ["CGGAGC", "TTAAGG"]


def test_get_codons_supports_reverse_and_reverse_complement_for_any_region():
    whitelists = {
        "S0": [
            {"codon": "GCCTCG"},
            {"codon": "AATTCC"},
        ],
    }
    reverse_only_region = Region(
        name="S0",
        index=0,
        codons=[],
        max_error_rate=0.0,
        indels=0,
        reverse=True,
    )
    reverse_complement_region = Region(
        name="S0",
        index=0,
        codons=[],
        max_error_rate=0.0,
        indels=0,
        reverse=True,
        complement=True,
    )

    assert get_codons(reverse_only_region, whitelists) == ["GCTCCG", "CCTTAA"]
    assert get_codons(reverse_complement_region, whitelists) == ["CGAGGC", "GGAATT"]


def test_generate_input_files_writes_transformed_region_adapters(tmp_path):
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
                "complement": True,
            },
            {
                "name": "S0",
                "max_error_rate": 0.0,
                "indels": 0,
                "reverse": True,
            },
            {
                "name": "C0",
                "max_error_rate": 0.0,
                "indels": 0,
                "reverse": True,
                "complement": True,
            },
        ],
        "whitelists": {
            "B0": [
                {"index": 0, "codon": "GCCTCG"},
                {"index": 1, "codon": "AATTCC"},
            ],
            "S0": [{"codon": "GCCTCG"}],
            "C0": [{"codon": "TTGGCC"}],
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
    c0_fastq = (
        tmp_path / "demo" / "demultiplex" / "cutadapt_input_files" / "2-C0.fastq"
    ).read_text()

    assert b0_fastq == ">0-B0.0\nCGGAGC\n>0-B0.1\nTTAAGG"
    assert s0_fastq == ">1-S0.0\nGCTCCG"
    assert c0_fastq == ">2-C0.0\nGGCCAA"
