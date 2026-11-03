# Automated Bulk RNA-seq Analysis Pipeline

An end-to-end, production-grade automated pipeline for Bulk RNA-seq processing and differential gene expression analysis. This project wraps industry-standard bioinformatics CLI tools inside a unified Python orchestrator and handles everything from raw sequencing files to interactive HTML diagnostics.

## Pipeline Architecture
The workflow is divided into modular steps, with outputs from each stage cleanly isolated:
1. **Quality Control (`FastQC`)** - Inspects raw sequencing reads for quality parameters.
2. **Trimming & Cleaning (`fastp`)** - Fast multi-threaded adapter clipping and base quality filtering.
3. **Splice-Aware Alignment (`HISAT2`)** - Maps clean paired-end reads to a reference genome.
4. **Quantification (`featureCounts`)** - Assigns mapped reads to genomic features (exons/genes).
5. **Differential Expression (`PyDESeq2`)** - Native Python implementation of the DESeq2 statistical model.
6. **Dashboard Aggregation (`MultiQC`)** - Compiles all metrics into a single interactive HTML document.

---

## Directory Structure
```text
├── setup.sh                     # Automates shell utility installations
├── requirements.txt             # Python packages layer
├── generate_mock_data.py        # Utility to build test sequences locally
├── pipeline.py         # Main execution pipeline script
├── fastq/                       # Input: Put your raw sample _R1/_R2 .fastq files here
├── reference/                   # Input: Reference genome.fa and annotations.gtf
└── output/             # Generated dynamically by the pipeline script
    ├── 1_quality_control/
    ├── 2_trimmed/
    ├── 3_aligned/
    ├── 4_quantification/
    ├── 5_differential_expression/ (Contains volcano_plot.png and de_results.csv)
    ├── dashboard/                 (Contains multiqc_report.html)
    └── logs/                      (Aggregated tool diagnostics)
```
### ⚠️ Permissions & Workspace Initialization
* **Write Access:** Ensure you have full write access to the project directory after cloning the repository.
* **Automated Setup:** When you run `complete_pipeline.py`, the orchestrator will automatically create all the required output directory structures dynamically. No manual folder configuration is needed.

---

## Getting Started

### 1. Prerequisites
Ensure you have [Miniconda](https://conda.io) or Anaconda installed on your Linux/macOS system.

### 2. Environment Setup
Run the automated initialization script to configure Conda channels, pull biological dependencies, and install Python libraries:
```bash
chmod +x setup.sh
./setup.sh
```

### 3. Generate Mock Test Data (Optional)
To verify everything is configured correctly without downloading hundreds of gigabytes of real genomic data, run the mock generator to create small valid sequence structures:
```bash
python generate_mock_data.py
```

### 4. Running the Pipeline
Execute the full workflow script. The pipeline automatically evaluates requirements, builds index baselines if missing, maps dependencies, runs statistics, and bundles logs:
```bash
python pipeline.py
```

---

## Outputs & Visualization
Upon successful completion, check the `output/` directory:
* **`output/5_differential_expression/de_results.csv`**: Complete spreadsheet table containing Log2 Fold Changes, p-values, and adjusted p-values per gene.
* **`output/5_differential_expression/volcano_plot.png`**: Volcano plot visualizing up-regulated vs down-regulated genes.
* **`output/dashboard/multiqc_report.html`**: A production-grade standalone web dashboard summarizing trimming yields, quality stats, and genomic mapping rates. Open this file directly in any web browser.
