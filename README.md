# Package `delt-core`
Core functionalities to work with DECL libraries.


## Software requirements and installation

Download Miniconda (select the installer that matches your operating system): https://docs.anaconda.com/miniconda. For macOS users: Choose the `pkg` installer, if you prefer a standard graphical installation process. If you are comfortable using the terminal, you can also use the `bash` installer. After installation, open the terminal and verify the installation by typing `conda --version`.

Use the following command to create a new environment called *del*: `conda create -n del` (type `y` and press **Enter** to proceed). You can list all your environments by typing `conda env list`. Activate the environment by typing `conda activate del`. The `delt-core` package will be installed only within the currently activated environment. **This also means that whenever you want to work with the** `delt-core` **package, you need to activate the environment first.** When an environment is activated, the name of the environment usually appears in parentheses at the beginning of the command prompt.

Generate a new SSH key by typing `ssh-keygen -t ed25519 -C "your_email@example.com"` in the terminal (enter your email address). When you're prompted to "Enter file in which to save the key", you can press **Enter** to accept the default file location. Start the ssh-agent by typing `eval "$(ssh-agent -s)"` and add your SSH private key to the ssh-agent by typing `ssh-add ~/.ssh/id_ed25519`. Copy the SSH public key to your clipboard by typing `pbcopy < ~/.ssh/id_ed25519.pub`. Go to GitHub, click on your profile picture (in the upper right corner), then click **Settings**. Go to **SSH and GPG keys** (in the sidebar) and click **New SSH key**. Add a title, paste the key in the corresponding field, and click **Add SSH key**.

Go to the repository on GitHub: https://github.com/DELTechnology/delt-core. Click the **Code** button, select **SSH**, and copy the URL to the clipboard. In the terminal, navigate to the directory where you want to clone the repository and type `git clone URL` (replace URL with the link in the clipboard). Type `cd delt-core` to navigate to the root folder of the package and run `conda install pip` and `pip install .` to install the package. Verify the installation by typing `delt-cli --help` (you should see a list of commands). Download the development version of Cutadapt: `pip install git+https://github.com/marcelm/cutadapt.git` (this command can be ignored once Cutadapt 4.10 is released).


## Installation (only for development purposes)

Navigate to the root folder of the package and run the following command:

```bash
pip install -e ".[dev,test]"
```

Create a `.env` file to store environment configurations, keys and secrets.
```bash
touch .env
```
These configurations can be accessed using the `python-dotenv` package.


## Example workflow

Initialize the folder structure and move the library, selection, and FASTQ files to the corresponding directories:
```bash
delt-cli init
mv /path/to/input.fastq.gz fastq_files
mv /path/to/library.xlsx libraries
mv /path/to/selection.xlsx selections
```

The folder structure should now be organized as follows:
```bash
.
├── fastq_files
│   └── input.fastq.gz
├── libraries
│   └── library.xlsx
└── selections
    └── selection.xlsx
```

Compute the SMILES and some chemical properties of a library:
```bash
delt-cli compute smiles libraries/library.xlsx
delt-cli compute properties libraries/smiles/library_smiles.txt.gz
delt-cli compute plot libraries/properties/properties_L1.txt.gz
```

Initialize the configuration file and demultiplex the FASTQ file (adjust the configuration file manually if needed):
```bash
delt-cli demultiplex init -f fastq_files/input.fastq.gz -l libraries/library.xlsx -s selections/selection.xlsx
delt-cli demultiplex run experiments/default-*/config.yml
```

Report and plot the results of the demultiplexing:
```bash
delt-cli qc report experiments/default-*
delt-cli qc plot experiments/default-*
```

Compare a set of target selections (e.g., ID 1-3) to a set of control selections (e.g., ID 4-6):
```bash
delt-cli normalize run experiments/default-*/config.yml '1 2 3' '4 5 6'
```


## Initialization

Initialize folder structure:
```bash
delt-cli init
```


## Computation

Compute SMILES of a library:
```bash
delt-cli compute smiles library1.xlsx
```

Compute SMILES of a hybridized library (the order of the libraries must match the final sequence in the 5'-to-3' direction, see templates/README.txt):
```bash
delt-cli compute smiles library1.xlsx library2.xlsx
```

Merge two library files to one Excel file (required for the demultiplexing of a hybridized library):
```bash
delt-cli compute merge library1.xlsx library2.xlsx
```

Compute chemical properties of a library:
```bash
delt-cli compute properties smiles/library1_smiles.txt.gz
```

Plot chemical properties of a library:
```bash
delt-cli compute plot properties/properties_L1.txt.gz
```


## Demultiplexing

Initialize folder structure for demultiplexing:
```bash
delt-cli demultiplex init
```

Run demultiplexing:
```bash
delt-cli demultiplex run experiments/default-*/config.yml
```

Create codon lists for the library specified in the configuration file:
```bash
delt-cli demultiplex create-lists experiments/default-*/config.yml
```

Create input files for Cutadapt:
```bash
delt-cli demultiplex create-cutadapt-input experiments/default-*/config.yml
```

Run bash script:
```bash
bash experiments/default-*/cutadapt_input_files/demultiplex.sh
```

Compute count tables for the final reads:
```bash
delt-cli demultiplex compute-counts experiments/default-*/config.yml experiments/default-*/cutadapt_output_files/reads_with_adapters.gz output_dir
```

Convert old structure file to new configuration file:
```bash
delt-cli demultiplex convert structure.txt
```


## Quality control

Plot codon hits:
```bash
delt-cli qc plot experiments/default-*
```

Print report:
```bash
delt-cli qc report experiments/default-*
```


## Normalization

Compare a set of target selections (e.g., ID 1-3) to a set of control selections (e.g., ID 4-6):
```bash
delt-cli normalize run experiments/default-*/config.yml '1 2 3' '4 5 6'
```


## Simulation

Create configuration file for simulation:
```bash
delt-cli simulate init
```

Generate reads (with or without erros):
```bash
delt-cli simulate run experiments/default-*/config.yml
```

The default initialization of the simulation generates a new library and a selection template that contain random codons. Alternatively, one can pass existing files using the following commands:
```bash
delt-cli simulate init -l libraries/library.xlsx -s selections/selection.xlsx -f fastq_files/input.fastq.gz -o fastq_files/input.fastq.gz
delt-cli simulate run experiments/default-*/config.yml
```
