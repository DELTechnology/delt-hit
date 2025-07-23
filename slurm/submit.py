import subprocess
import textwrap

def submit_job(args, debug: bool = False):
    script = f"""#!/bin/bash
#SBATCH --job-name=delt
#SBATCH --output=/work/FAC/FBM/DBC/mrapsoma/prometex/logs/adrianom/delt-%j.log
#SBATCH --error=/work/FAC/FBM/DBC/mrapsoma/prometex/logs/adrianom/delt-%j.err

#SBATCH --account=mrapsoma_prometex

#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# Load environment
source ~/envs/delt-core/bin/activate

cd /users/amarti51/projects/delt-core || exit

python {args}
"""
    script = textwrap.dedent(script)

    if debug:
        print(script)
        return

    proc = subprocess.run(
        ["sbatch"],
        input=script,
        text=True,
        check=True,
        capture_output=True,
    )
    job_id = proc.stdout.strip()
    print(f"Submitted array job {job_id}")

debug = False
args = '/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-core/scripts/01-embed-ChEMBL.py'
submit_job(args=args, debug=debug)
