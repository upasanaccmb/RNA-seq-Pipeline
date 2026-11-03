#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "=========================================================="
echo " Starting Bulk RNA-Seq Pipeline Environment Setup         "
echo "=========================================================="

# 1. Verify conda is installed
if ! command -v conda &> /dev/null; then
    echo "[ERROR] Conda could not be found. Please install Miniconda or Anaconda first."
    exit 1
fi

echo "[-] Configuring Conda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

echo "[-] Installing CLI Bioinformatics dependencies via Conda..."
# FastQC: Quality Control
# fastp: Trimming & Cleaning
# hisat2: Splice-aware alignment
# subread: Includes featureCounts for quantification
conda install -y fastqc fastp hisat2 subread

echo "[-] Installing Python requirements via pip..."
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo "[WARNING] requirements.txt not found! Skipping pip installation."
fi

echo "=========================================================="
echo "[SUCCESS] Environment is fully configured and ready!"
echo "=========================================================="
