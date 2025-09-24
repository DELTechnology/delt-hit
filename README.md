# 🧬 `delt-core`
> Core functionalities to work with DNA-encoded chemical libraries.

## 🚀 Installation

This guide provides instructions for setting up `delt-core` for both regular users and developers.

### Prerequisites

Before you begin, make sure you have the following installed:

#### 1. Conda
We recommend using the [Miniconda](https://docs.anaconda.com/miniconda) package manager to create an isolated environment for this project. This ensures that all dependencies are managed correctly.
- [Download and install Miniconda](https://docs.anaconda.com/miniconda#latest-miniconda-installer-links) for your operating system.
- After installation, you should be able to use the `conda` command in your terminal.

#### 2. R Environment
Some analysis features in `delt-core` (like enrichment analysis with `edgeR`) depend on R.
- **Install R:** Download and install R from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).
- **Install R Packages:** Once R is installed, open an R console and run the following commands to install the required packages:
    ```R
    # Install tidyverse and GGally from CRAN
    install.packages(c("tidyverse", "GGally"))

    # Install BiocManager
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    # Install edgeR and limma from Bioconductor
    BiocManager::install(c("edgeR", "limma"))
    ```

### 🧑‍🔬 User Installation

This is the recommended way for most users.

1.  **Create and activate a Conda environment:**
    ```bash
    conda create -n del python=3.11 -y
    conda activate del
    ```
    > 💡 Always activate this environment (`conda activate del`) before using `delt-core`.

2.  **Install `delt-core`:**
    Install the package directly from GitHub using `pip`:
    ```bash
    pip install git+https://github.com/DELTechnology/delt-core.git
    ```
    > **Note:** The `delt-core` package is under active development. To get the latest version of `cutadapt` required by this package, please run `pip install git+https://github.com/marcelm/cutadapt.git` (this command can be ignored once Cutadapt 4.10 is released).

3.  **Verify Installation:**
    Check that the CLI is working:
    ```bash
    delt-cli --help
    ```
    You should see a list of available commands.

### 👩‍💻 Developer Installation

If you want to contribute to the development of `delt-core`, follow these steps.

1.  **Configure SSH for GitHub:**
    Make sure you have an [SSH key added to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) to clone the repository.

2.  **Clone the Repository:**
    ```bash
    git clone git@github.com:DELTechnology/delt-core.git
    cd delt-core
    ```

3.  **Create and activate the Conda environment:**
    ```bash
    conda create -n del-dev python=3.11 -y
    conda activate del-dev
    ```

4.  **Install in Editable Mode:**
    Install the package with all development and testing dependencies:
    ```bash
    pip install -e ".[dev,test]"
    ```
    > 🔧 This "editable" install means that any changes you make to the source code will be immediately reflected when you run the `delt-cli` command.

5.  **(Optional) Install `pigz` for parallel processing:**
    For faster demultiplexing on macOS, install `pigz` using [Homebrew](https://brew.sh/):
    ```bash
    brew install pigz
    ```

## 🧪 Example Workflow

Here is a typical workflow for using `delt-core`:

1.  **Initialize Configuration:**
    Create a `config.yaml` file from an Excel library file. This file defines the experiment, selections, and library information.
    ```bash
    delt-cli demultiplex init --excel_path /path/to/library.xlsx
    ```

2.  **Run Demultiplexing:**
    Run the entire demultiplexing pipeline based on your configuration.
    ```bash
    delt-cli demultiplex run --config_path /path/to/config.yaml
    ```

3.  **Define Analysis Groups:**
    After demultiplexing, group selections for analysis.
    ```bash
    delt-cli analyse add \
    --config_path /path/to/config.yaml \
    --name=test-1 \
    --selections='["SEL1", "SEL2", "SEL3"]'
    ```

4.  **Calculate Enrichment:**
    Calculate enrichment for the defined groups using different methods.
    ```bash
    delt-cli analyse enrichment --config_path /path/to/config.yaml --name=test-1 --method=counts
    delt-cli analyse enrichment --config_path /path/to/config.yaml --name=test-1 --method=edgeR
    ```

5.  **Work with the Library:**
    Enumerate the library, compute properties, and generate representations.
    ```bash
    # Enumerate all molecules in the library
    delt-cli library enumerate --config_path /path/to/config.yaml

    # Compute chemical properties
    delt-cli library properties --config_path /path/to/config.yaml

    # Generate molecular fingerprints
    delt-cli library represent --method=morgan --config_path /path/to/config.yaml
    ```

6.  **Launch Dashboard:**
    Explore the results interactively in a web-based dashboard.
    ```bash
    delt-cli dashboard \
    --config_path /path/to/config.yaml \
    --counts_path /path/to/selections/SELECTION_NAME/counts.txt
    ```

## 💻 CLI Reference

### 📚 Library
Commands for library enumeration, and chemical property and representation calculation.

- **`enumerate`**: Generates the full library of molecules from the reaction steps defined in the configuration file.
  ```bash
  delt-cli library enumerate --config_path <path/to/config.yaml>
  ```
- **`properties`**: Calculates a set of chemical properties for the enumerated library.
  ```bash
  delt-cli library properties --config_path <path/to/config.yaml>
  ```
- **`represent`**: Generates molecular representations (fingerprints) for the library.
  ```bash
  delt-cli library represent --config_path <path/to/config.yaml> --method <METHOD>
  ```
  - `<METHOD>` can be `morgan` or `bert`.

### ✂️ Demultiplexing
Commands for demultiplexing FASTQ files and obtaining read counts.

- **`init`**: Creates a `config.yaml` file from a library Excel file.
  ```bash
  delt-cli demultiplex init --excel_path <path/to/library.xlsx>
  ```
- **`run`**: Runs the entire demultiplexing workflow, including running Cutadapt and computing counts.
  ```bash
  delt-cli demultiplex run --config_path <path/to/config.yaml>
  ```
- Other steps (`prepare`, `process`, `report`, `qc`) can be run individually for debugging or custom workflows.

### 📊 Analysis
Commands for analyzing demultiplexed data, such as performing enrichment analysis.

- **`add`**: Defines a named group of selections for downstream analysis.
  ```bash
  delt-cli analyse add --config_path <path/to/config.yaml> --name <group_name> --selections '["SEL1", "SEL2", ...]' 
  ```
- **`enrichment`**: Performs enrichment analysis on a defined analysis group.
  ```bash
  delt-cli analyse enrichment --config_path <path/to/config.yaml> --name <group_name> --method <METHOD>
  ```
  - `<METHOD>` can be `counts`, `edgeR`, or `DESeq2`.

### 📈 Dashboard
Launch an interactive dashboard for data visualization.

- **`dashboard`**: Starts a web-based dashboard to interactively explore counts data for a given selection.
  ```bash
  delt-cli dashboard --config_path <path/to/config.yaml> --counts_path <path/to/counts.txt>
  ```