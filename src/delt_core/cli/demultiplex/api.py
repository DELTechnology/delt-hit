import json
import subprocess
from pathlib import Path

from delt_core.demultiplex.parser import config_from_excel
from delt_core.demultiplex.postprocess import get_counts, save_counts
from delt_core.demultiplex.preprocess import generate_input_files
from delt_core.utils import write_yaml, read_yaml
from loguru import logger

class Demultiplex:

    def prepare(self, *, config_path: Path, fast_dev_run: bool = False):
        exec_path = generate_input_files(config_path=config_path, fast_dev_run=fast_dev_run)
        logger.info(f"Executable created at {exec_path}")

    def process(self, *, config_path: Path):
        config = read_yaml(config_path)
        save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
        name = config['experiment']['name']

        output_dir = save_dir / name / 'demultiplex' / 'cutadapt_output_files'
        input_path = output_dir / 'reads_with_adapters.gz'

        num_reads = json.load(open(
            sorted(output_dir.glob('*.cutadapt.json'))[-1]
        ))['read_counts']['output']

        counts = get_counts(input_path=input_path, num_reads=num_reads)

        ids_to_name = {tuple(item['ids']): k for k, item in config['selections'].items()}
        output_dir = save_dir / name / 'selections'
        save_counts(counts, output_dir=output_dir, ids_to_name=ids_to_name)

    def report(self, *, config_path: Path):
        from delt_core.quality_control.report import print_report
        config = read_yaml(config_path)
        save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
        name = config['experiment']['name']

        output_dir = save_dir / name / 'demultiplex' / 'cutadapt_output_files'
        save_path = save_dir / name / 'qc' / 'report.txt'
        save_path.parent.mkdir(parents=True, exist_ok=True)
        print_report(output_dir=output_dir, save_path=save_path)

    def qc(self, *, config_path: Path):
        from delt_core.quality_control.plot_codon_hits import plot_hits

        config = read_yaml(config_path)
        save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
        name = config['experiment']['name']

        output_dir = save_dir / name / 'demultiplex' / 'cutadapt_output_files'
        save_dir = save_dir / name / 'qc'
        save_dir.mkdir(parents=True, exist_ok=True)
        plot_hits(output_dir=output_dir, save_dir=save_dir)

    def run(self, *, config_path: Path, fast_dev_run: bool = False):
        exec_path = generate_input_files(config_path=config_path, fast_dev_run=fast_dev_run)
        subprocess.run(['bash', exec_path])