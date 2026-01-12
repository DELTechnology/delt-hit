from delt_core.quality_control.report import print_report
from pathlib import Path

output_dir = Path('/Volumes/T7/experiments/experiment-GB-minimal/demultiplex/cutadapt_output_files')
print_report(output_dir=output_dir, save_path=output_dir / 'report.txt')