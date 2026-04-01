"""
DELT-Hit Snakemake pipeline
============================
Chains all delt-hit CLI commands from raw FASTQ to enrichment results.

Usage:
    snakemake --cores 4                  # run full pipeline
    snakemake --cores 4 rule_name        # run a specific rule
    snakemake -n                         # dry-run (show what would execute)
    snakemake -n --quiet                 # dry-run, concise output
    snakemake --dag | dot -Tpng > dag.png  # visualise the dependency graph
    snakemake --forceall --cores 4       # re-run everything regardless of timestamps

Configuration:
    This pipeline reads two YAML files:
      - config.yaml       : demultiplex config (experiment, selections, structure, whitelists)
      - analysis.yaml     : enrichment config  (experiments with counts_path per selection)
    Both are read by the corresponding delt-hit commands. Snakemake reads config.yaml
    to resolve file paths; analysis.yaml is passed directly to delt-hit analyse enrichment.

Optional naive library analysis:
    Add `naive_fastq_path: /path/to/naive.fastq.gz` under `experiment:` in config.yaml
    to enable the naive demultiplex and library assess rules.
"""

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
configfile: "config.yaml"

# Convenience aliases derived from the delt-hit config schema
SAVE_DIR   = config["experiment"]["save_dir"]
NAME       = config["experiment"]["name"]
EXP_DIR    = f"{SAVE_DIR}/{NAME}"
FASTQ      = config["experiment"]["fastq_path"]
NUM_CORES  = config["experiment"].get("num_cores", 1)
SELECTIONS = list(config["selections"].keys())

# Optional: separate analysis config for enrichment step
ANALYSIS_CONFIG = config.get("analysis_config", "analysis.yaml")

# Optional naive library FASTQ (enables rule naive_demultiplex + rule assess)
NAIVE_FASTQ = config["experiment"].get("naive_fastq_path", None)

# ---------------------------------------------------------------------------
# Rule all — default target; collecting the terminal outputs of every branch
# ---------------------------------------------------------------------------
# Branches:
#   A. Demultiplex:  FASTQ → reads_with_adapters.gz → counts.txt per selection
#   B. QC:          cutadapt outputs → qc/report.txt
#   C. Library:     config → library.parquet → properties.parquet + morgan.npz
#   D. Enrichment:  counts → R script (Rscript must be run separately)
# Optional:
#   E. Naive:       naive FASTQ → naive counts + library assess report

_required_targets = [
    # Demultiplex branch
    f"{EXP_DIR}/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
    expand(f"{EXP_DIR}/selections/{{sel}}/counts.txt", sel=SELECTIONS),
    # QC branch
    f"{EXP_DIR}/qc/report.txt",
    # Library branch
    f"{EXP_DIR}/library/library.parquet",
    f"{EXP_DIR}/library/properties/properties.parquet",
    f"{EXP_DIR}/representations/morgan.npz",
    # Enrichment branch (R script is the final Python-produced artefact)
    f"{EXP_DIR}/counts/enrichment_counts.R",
]

_optional_targets = (
    [
        f"{EXP_DIR}/naive/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
        f"{EXP_DIR}/naive/qc/report.txt",
    ]
    if NAIVE_FASTQ
    else []
)

rule all:
    """
    Default target. Runs the complete pipeline.
    All branches are independent of each other after `init` and can be
    executed in parallel with --cores N.
    """
    input:
        _required_targets + _optional_targets


# ---------------------------------------------------------------------------
# Rule init — Excel workbook → YAML config
# ---------------------------------------------------------------------------
# Input:  Excel file (templates/library.xlsx or custom)
# Output: {EXP_DIR}/config.yaml
# Note:   Only needed when starting from an Excel spec. If you already have
#         a config.yaml, skip this rule — Snakemake will not re-run it.
# ---------------------------------------------------------------------------
rule init:
    """Convert an Excel workbook to a delt-hit YAML config file."""
    input:
        excel = config.get("excel_path", "templates/library.xlsx"),
    output:
        cfg = f"{EXP_DIR}/config.yaml",
    log:
        f"{EXP_DIR}/logs/init.log",
    shell:
        "delt-hit init --excel_path={input.excel} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Rule demultiplex — FASTQ → demultiplexed reads + counts per selection
# ---------------------------------------------------------------------------
# This rule runs the full demultiplex pipeline:
#   1. Generates cutadapt input files and demultiplex.sh  (prepare)
#   2. Executes demultiplex.sh via bash                   (run)
#   3. Counts reads per selection and writes counts.txt    (process)
#
# Parallelism: cutadapt itself is multi-threaded (num_cores from config).
#              Snakemake threads= here matches that so the scheduler can
#              reserve the right number of cores.
#
# Input:  FASTQ file + config.yaml
# Output: reads_with_adapters.gz  (sentinel for demultiplex run)
#         selections/{name}/counts.txt  (one per selection)
# ---------------------------------------------------------------------------
rule demultiplex:
    """Run cutadapt demultiplexing and produce per-selection count tables."""
    input:
        fastq  = FASTQ,
        config = "config.yaml",
    output:
        reads  = f"{EXP_DIR}/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
        counts = expand(f"{EXP_DIR}/selections/{{sel}}/counts.txt", sel=SELECTIONS),
    log:
        f"{EXP_DIR}/logs/demultiplex.log",
    threads:
        NUM_CORES
    shell:
        """
        delt-hit demultiplex run \
            --config_path={input.config} \
            > {log} 2>&1

        delt-hit demultiplex process \
            --config_path={input.config} \
            --as_files=false \
            --sort_by_counts=true \
            >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule qc — cutadapt outputs → QC report + codon-hit plots
# ---------------------------------------------------------------------------
# Runs after demultiplex (needs cutadapt JSON files).
# Parallel with: rule enumerate, rule enrichment.
#
# Input:  reads_with_adapters.gz (sentinel — actual inputs are the JSON logs
#         in the same directory, which cutadapt wrote during demultiplex)
# Output: qc/report.txt
# ---------------------------------------------------------------------------
rule qc:
    """Generate the cutadapt summary report and codon-hit QC plots."""
    input:
        reads  = f"{EXP_DIR}/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
        config = "config.yaml",
    output:
        report = f"{EXP_DIR}/qc/report.txt",
    log:
        f"{EXP_DIR}/logs/qc.log",
    shell:
        """
        delt-hit demultiplex report \
            --config_path={input.config} \
            > {log} 2>&1

        delt-hit demultiplex qc \
            --config_path={input.config} \
            >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule enumerate — config → combinatorial library (SMILES)
# ---------------------------------------------------------------------------
# Independent of the demultiplex branch — can run in parallel with it.
#
# Input:  config.yaml  (reads library/catalog/whitelists sections)
# Output: library/library.parquet
# ---------------------------------------------------------------------------
rule enumerate:
    """Enumerate all BB combinations and compute product SMILES via RDKit."""
    input:
        config = "config.yaml",
    output:
        library = f"{EXP_DIR}/library/library.parquet",
    log:
        f"{EXP_DIR}/logs/enumerate.log",
    shell:
        """
        delt-hit library enumerate \
            --config_path={input.config} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule properties — library.parquet → molecular property distributions
# ---------------------------------------------------------------------------
# Parallel with: rule represent.
#
# Input:  library/library.parquet
# Output: library/properties/properties.parquet
# ---------------------------------------------------------------------------
rule properties:
    """Compute RDKit molecular properties (MW, logP, QED, …) for the library."""
    input:
        library = f"{EXP_DIR}/library/library.parquet",
        config  = "config.yaml",
    output:
        props = f"{EXP_DIR}/library/properties/properties.parquet",
    log:
        f"{EXP_DIR}/logs/properties.log",
    threads:
        max(1, NUM_CORES // 2)   # property calc is embarrassingly parallel; share cores
    shell:
        """
        delt-hit library properties \
            --config_path={input.config} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule represent — library.parquet → Morgan fingerprints
# ---------------------------------------------------------------------------
# Parallel with: rule properties.
# Use --method=bert for transformer embeddings (requires [ml] extras + GPU).
#
# Input:  library/library.parquet
# Output: representations/morgan.npz
# ---------------------------------------------------------------------------
rule represent:
    """Compute Morgan fingerprints for the library (sparse .npz format)."""
    input:
        library = f"{EXP_DIR}/library/library.parquet",
        config  = "config.yaml",
    output:
        fps = f"{EXP_DIR}/representations/morgan.npz",
    log:
        f"{EXP_DIR}/logs/represent.log",
    shell:
        """
        delt-hit library represent \
            --config_path={input.config} \
            --method=morgan \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule enrichment — counts.txt per selection → enrichment R script
# ---------------------------------------------------------------------------
# NOTE: This rule generates an R script but does NOT execute it.
#       After this rule completes, run:
#           Rscript {EXP_DIR}/counts/enrichment_counts.R
#       or:
#           Rscript {EXP_DIR}/edgeR/enrichment_edgeR.R
#
# The enrichment command reads from analysis.yaml (different schema from
# config.yaml). analysis.yaml lists experiments with their selections and
# the path to each selection's counts.txt.
#
# Input:  counts.txt per selection + analysis.yaml
# Output: enrichment_counts.R  (the generated R script)
# ---------------------------------------------------------------------------
rule enrichment:
    """Generate enrichment analysis R script (edgeR or simple counts method)."""
    input:
        counts = expand(f"{EXP_DIR}/selections/{{sel}}/counts.txt", sel=SELECTIONS),
        config = ANALYSIS_CONFIG,
    output:
        rscript = f"{EXP_DIR}/counts/enrichment_counts.R",
    log:
        f"{EXP_DIR}/logs/enrichment.log",
    params:
        name   = NAME,
        method = config.get("enrichment_method", "counts"),
    shell:
        """
        delt-hit analyse enrichment \
            --config_path={input.config} \
            --name={params.name} \
            --method={params.method} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Optional: Naive library demultiplex
# ---------------------------------------------------------------------------
# Activated when experiment.naive_fastq_path is set in config.yaml.
# Produces a separate counts table for the unselected (naive) library,
# which feeds into synthesis QC and normalization of selection counts.
#
# Input:  naive FASTQ + config.yaml
# Output: naive/demultiplex/cutadapt_output_files/reads_with_adapters.gz
#         naive/qc/report.txt
# ---------------------------------------------------------------------------
if NAIVE_FASTQ:
    rule naive_demultiplex:
        """
        Demultiplex the naive (unselected) library FASTQ.
        Naive counts reveal synthesis yield per building block and allow
        normalization of selection enrichment scores.
        """
        input:
            fastq  = NAIVE_FASTQ,
            config = "config.yaml",
        output:
            reads = f"{EXP_DIR}/naive/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
        log:
            f"{EXP_DIR}/logs/naive_demultiplex.log",
        threads:
            NUM_CORES
        shell:
            # Re-use the same demultiplex command pointed at a temporary config
            # with naive_fastq_path substituted for fastq_path.
            # A helper script writes the modified config; adjust if you have
            # a dedicated naive config file instead.
            """
            python - <<'EOF'
import yaml, pathlib, shutil

cfg = yaml.safe_load(open("{input.config}"))
cfg["experiment"]["fastq_path"] = "{input.fastq}"
cfg["experiment"]["save_dir"] = str(pathlib.Path("{EXP_DIR}") / "naive" / "..")
cfg["experiment"]["name"] = "{NAME}/naive"
pathlib.Path("{EXP_DIR}/naive").mkdir(parents=True, exist_ok=True)
yaml.dump(cfg, open("{EXP_DIR}/naive/config.yaml", "w"), sort_keys=False)
EOF
            delt-hit demultiplex run \
                --config_path={EXP_DIR}/naive/config.yaml \
                > {log} 2>&1
            """

    rule naive_qc:
        """QC report for the naive library demultiplex run."""
        input:
            reads  = f"{EXP_DIR}/naive/demultiplex/cutadapt_output_files/reads_with_adapters.gz",
            config = f"{EXP_DIR}/naive/config.yaml",
        output:
            report = f"{EXP_DIR}/naive/qc/report.txt",
        log:
            f"{EXP_DIR}/logs/naive_qc.log",
        shell:
            """
            delt-hit demultiplex report \
                --config_path={input.config} \
                > {log} 2>&1

            delt-hit demultiplex qc \
                --config_path={input.config} \
                >> {log} 2>&1
            """
