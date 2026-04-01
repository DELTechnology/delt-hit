# delt-hit Dockerfile
#
# NOTE: python:3.12-slim is used instead of 3.11-slim because the source code
# uses itertools.batched (added in Python 3.12).
#
# Build:  docker build -t delt-hit .
# Run:    docker run --rm -v $(pwd)/data:/workspace/data delt-hit --help

FROM python:3.12-slim

# ---------------------------------------------------------------------------
# System dependencies
# ---------------------------------------------------------------------------
# - gcc / g++ / make: required to compile C extensions (e.g. rdkit, cutadapt)
# - git: some uv/pip installs resolve git dependencies
# - r-base: R runtime needed for edgeR analysis scripts
# - libcurl4-openssl-dev / libssl-dev: required by R package installer
# Packages are installed in a single RUN layer and the apt cache is cleared
# immediately to keep the image small.
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        make \
        git \
        r-base \
        libcurl4-openssl-dev \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Install uv (fast Python package manager)
# ---------------------------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# ---------------------------------------------------------------------------
# Set working directory
# ---------------------------------------------------------------------------
WORKDIR /workspace

# ---------------------------------------------------------------------------
# Install Python dependencies (cached layer — only rebuilt when pyproject.toml
# or uv.lock change, not when source code changes)
# ---------------------------------------------------------------------------
COPY pyproject.toml ./
# Copy uv.lock only if it exists (use a wildcard so the build doesn't fail
# when the file is absent)
COPY uv.loc[k] ./

# Install all extras so the full CLI is available; no editable install yet
# (source code is copied in the next step)
RUN uv sync --all-extras --no-install-project

# ---------------------------------------------------------------------------
# Install R packages: edgeR from Bioconductor
# ---------------------------------------------------------------------------
# This layer is slow (~5 min) but cached as long as nothing above changes.
RUN Rscript -e '\
    if (!require("BiocManager", quietly = TRUE)) \
        install.packages("BiocManager", repos = "https://cloud.r-project.org"); \
    BiocManager::install("edgeR", ask = FALSE, update = FALSE)'

# ---------------------------------------------------------------------------
# Copy source code and install the package
# ---------------------------------------------------------------------------
COPY src/ ./src/

# Install the package itself (editable so the mounted /workspace/src works)
RUN uv sync --all-extras

# ---------------------------------------------------------------------------
# Default entrypoint
# ---------------------------------------------------------------------------
ENTRYPOINT ["delt-hit"]
