# Package `delt-core`
Core functionalities to work with DECL libraries

# Installation

Install the package for development purposes.
Navigate to the root folder of the package and run the following command:

```bash
pip install -e ".[dev,test]"
```

Create a `.env` file to store environment configurations, keys and secrets.
```bash
touch .env
```
These configurations can be accessed using the `python-dotenv` package.


# CLI

In a terminal, run the following commands to test the CLI:

```bash
delt-cli --help
delt-cli welcome # Welcome to DEL Technology!
delt-cli process # Processing...
```

# DECL Simulation

## FASTQ Simulations
The simulations uses the same input as the current evaulation program for the fastq files.
The simulation is configured by a `simulation/config.json` file that specifies the following keys:

- path_to_struct_file: Path to the struct file that specifies the structure of the fastq files
- number_of_reads: Number of reads to simulate
- errors: List of errors to simulate. Each element in the list is a dict of `type` that  defines the error type and a
    `kwargs` dict of arguments passed to the error function. The error types are:
  - all_barcode_regions_single_position: Simulates an error in all barcode regions at a single position

The simulated data is stored in the fastq file defined in the struct file.

```
python simulate-fastq.py <PATH_TO_SIMULATION_CONFIGURATION_FILE>
python simulate-fastq.py simulation/config.json
```

# TODO
- Illumina adapter trimmiing
- Check if after adapater removal the reads are aligned
- Create pipeline to sequentially demultiplex the fastq files based on barcode regions
