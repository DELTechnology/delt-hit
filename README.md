# 🧬 `DELT-Hit`

Welcome to `delt-hit`! An end-to-end computational framework for DNA-encoded chemical library analysis.

## 🚀 Installation

This guide provides instructions for setting up `delt-hit` for both regular users and developers.
Use either **[Pixi](https://pixi.sh/latest/)** (recommended) or **[Conda](https://docs.anaconda.com/miniconda)** to manage the environment.

### 🧑‍🔬 User Installation

Clone the repository:

```bash
git clone https://github.com/DELTechnology/delt-hit.git
cd delt-hit
```

#### Option A: Pixi (recommended)

```bash
pixi install
pixi run delt-hit --help
```

#### Option B: Conda

```bash
conda create -n delt-hit python=3.12 -y
conda activate delt-hit
conda install pygraphviz -y
pip install .
delt-hit --help
```

#### R dependencies (optional)

Required only for enrichment analysis with `edgeR`:

```R
install.packages(c("tidyverse", "GGally"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma"))
```

### 👩‍💻 Developer Installation

Same as the user installation above, but pass the `-e` flag to install in editable mode with dev dependencies:

- **Pixi:** `pixi install -e dev`
- **Conda:** `pip install -e ".[dev,test]"` instead of `pip install .`

## 🧪 Example Workflow

A complete end-to-end example is available in [`supporting_material/experiments/example-single-stranded/`](supporting_material/experiments/example-single-stranded/); see [`supporting_material/README.md`](supporting_material/README.md) for full instructions.

Here is a typical workflow for using `delt-hit`:

1. **Initialize Configuration:**
   Create a `config.yaml` file from an Excel library file. This file defines the experiment, selections, and library
   information.
   ```bash
   delt-hit init --excel_path /path/to/library.xlsx
   ```

2. **Run Demultiplexing:**
   Run the entire demultiplexing pipeline based on your configuration. This includes preparing scripts, running
   `cutadapt`, and processing the results.
   ```bash
   delt-hit demultiplex run --config_path /path/to/config.yaml
   ```

3. **Define Analysis Groups:**
   After demultiplexing, define analysis groups by editing your `config.yaml` file. Add an `analyses` section to group
   selections for comparison. For example:
   ```yaml
   experiments:
   - name: condition_vs_control
     save_dir: campaign/analysis
     selections:
       - name: AG24_4
         counts_path: campaign/selections/AG24_4/counts.txt
         group: control
       - name: AG24_5
         counts_path: campaign/selections/AG24_5/counts.txt
         group: control
       - name: AG24_6
         counts_path: campaign/selections/AG24_6/counts.txt
         group: control
       - name: AG24_13
         counts_path: campaign/selections/AG24_13/counts.txt
         group: condition
       - name: AG24_14
         counts_path: campaign/selections/AG24_14/counts.txt
         group: condition
       - name: AG24_15
         counts_path: campaign/selections/AG24_15/counts.txt
         group: condition
     ```

4. **Calculate Enrichment:**
   Calculate enrichment for the defined groups using different methods. The `--name` argument must correspond to a group
   you defined in your `config.yaml`.
   ```bash
   # Using simple counts
   delt-hit analyse enrichment --config_path /path/to/config.yaml --name=condition_vs_control --method=counts

   # Using edgeR for more sensitive statistical analysis
   delt-hit analyse enrichment --config_path /path/to/config.yaml --name=condition_vs_control --method=edgeR
   ```

5. **Work with the Library:**
   Enumerate the library, compute properties, and generate representations.
   ```bash
   # Enumerate all molecules in the library
   delt-hit library enumerate --config_path /path/to/config.yaml

   # Enumerate only the top observed combinations from a demultiplex counts file
   delt-hit library enumerate \
   --config_path /path/to/config.yaml \
   --counts_path /path/to/selections/SELECTION_NAME/counts.txt \
   --top_n 1000 \
   --library_name observed_hits

   # Compute chemical properties
   delt-hit library properties --config_path /path/to/config.yaml

   # Compute chemical properties for a named filtered library
   delt-hit library properties \
   --config_path /path/to/config.yaml \
   --library_name observed_hits

   # Generate molecular fingerprints (e.g., Morgan)
   delt-hit library represent --method=morgan --config_path /path/to/config.yaml
   ```

6. **Launch Dashboard:**
   Explore the results interactively in a web-based dashboard.
   ```bash
   delt-hit dashboard \
   --config_path /path/to/config.yaml \
   --counts_path /path/to/selections/SELECTION_NAME/counts.txt
   ```

## 📚 Documentation

For a codebase overview and a detailed CLI reference, see:

- [Codebase overview](documentation/overview.md)
- [CLI guide](documentation/cli.md)

The original protocol description lives in `protocols.pdf`.

## 💻 CLI Reference

For the most up-to-date CLI details and output locations, use the [CLI guide](documentation/cli.md).

### `init`

Initializes a project by creating a `config.yaml` from a standardized Excel file.

```bash
delt-hit init --excel_path <path/to/library.xlsx>
```

### `library`

Commands for library enumeration, and chemical property and representation calculation.

- **`enumerate`**: Generates the full library of molecules from the reaction steps defined in the configuration file.
  ```bash
  delt-hit library enumerate --config_path <path/to/config.yaml>
  ```
  For troubleshooting, `--debug` can write per-combination reaction-graph PNGs into the library output directory. Use `--debug invalid --errors ignore` to keep enumerating while capturing combinations that fail during reaction execution.
  ```bash
  delt-hit library enumerate \
  --config_path <path/to/config.yaml> \
  --debug invalid \
  --errors ignore
  ```
  You can also enumerate only the top observed barcode combinations from a demultiplex counts file:
  ```bash
  delt-hit library enumerate \
  --config_path <path/to/config.yaml> \
  --counts_path <path/to/selections/SELECTION_NAME/counts.txt> \
  --top_n 1000 \
  --library_name observed_hits
  ```
- **`properties`**: Calculates a set of chemical properties for the enumerated library.
  ```bash
  delt-hit library properties --config_path <path/to/config.yaml>
  ```
  You can also compute properties for a named filtered library:
  ```bash
  delt-hit library properties \
  --config_path <path/to/config.yaml> \
  --library_name observed_hits
  ```
- **`represent`**: Generates molecular representations (fingerprints) for the library.
  ```bash
  delt-hit library represent --config_path <path/to/config.yaml> --method <METHOD>
  ```
    - `<METHOD>` can be `morgan` or `bert`.

### `demultiplex`

Commands for demultiplexing FASTQ files and obtaining read counts.

- **`run`**: Runs the entire demultiplexing workflow, including running Cutadapt and computing counts.
  ```bash
  delt-hit demultiplex run --config_path <path/to/config.yaml>
  ```
- **`prepare`**: Prepares the `cutadapt` input files and executable script without running them.
  ```bash
  delt-hit demultiplex prepare --config_path <path/to/config.yaml>
  ```
- **`process`**: Computes counts from the output of a `cutadapt` run.
  ```bash
  delt-hit demultiplex process --config_path <path/to/config.yaml>
  ```
- **`report`**: Generates a text report summarizing demultiplexing statistics.
  ```bash
  delt-hit demultiplex report --config_path <path/to/config.yaml>
  ```
- **`qc`**: Generates quality control plots from the demultiplexing results.
  ```bash
  delt-hit demultiplex qc --config_path <path/to/config.yaml>
  ```

### `analyse`

Commands for analyzing demultiplexed data, such as performing enrichment analysis.

- **`enrichment`**: Performs enrichment analysis on an analysis group defined in the configuration file.
  ```bash
  delt-hit analyse enrichment --config_path <path/to/config.yaml> --name <group_name> --method <METHOD>
  ```
    - Analysis groups must be defined under the `analyses` key in your `config.yaml`.
    - `<group_name>` refers to a key under the `analyses` section.
    - `<METHOD>` can be `counts` or `edgeR`.

### `dashboard`

Launches an interactive dashboard for data visualization.

- **`dashboard`**: Starts a web-based dashboard to interactively explore counts data for a given selection.
  ```bash
  delt-hit dashboard --config_path <path/to/config.yaml> --counts_path <path/to/counts.txt>
  ```
