import os
import subprocess
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# =====================================================================
# CONFIGURATION & METADATA ARCHIVE
# =====================================================================
FASTQ_DIR = "./fastq"
REF_DIR = "./reference"
GENOME_FA = f"{REF_DIR}/genome.fa"
ANNOTATION_GTF = f"{REF_DIR}/annotations.gtf"

# Output Pipeline Directory Structure
BASE_OUT = "./output"
QC_DIR = f"{BASE_OUT}/1_quality_control"
TRIM_DIR = f"{BASE_OUT}/2_trimmed"
ALIGN_DIR = f"{BASE_OUT}/3_aligned"
QUANT_DIR = f"{BASE_OUT}/4_quantification"
DE_DIR = f"{BASE_OUT}/5_differential_expression"
LOGS_DIR = f"{BASE_OUT}/logs"          # Central repository for MultiQC parsing
DASHBOARD_DIR = f"{BASE_OUT}/dashboard" # Final HTML location

# Explicit Metadata Sample Sheet (Experimental Design)
SAMPLE_INFO = pd.DataFrame({
    "sample_id": ["sample1", "sample2", "sample3", "sample4"],
    "condition": ["Control", "Control", "Treatment", "Treatment"]
}).set_index("sample_id")


def initialize_workspace():
    """Ensures clean infrastructure output paths."""
    for folder in [QC_DIR, TRIM_DIR, ALIGN_DIR, QUANT_DIR, DE_DIR, LOGS_DIR, DASHBOARD_DIR]:
        os.makedirs(folder, exist_ok=True)
    print("[+] Pipeline directory infrastructure initialized.")


# =====================================================================
# STEP 1: QUALITY CONTROL (FastQC)
# =====================================================================
def run_quality_control():
    print("\n[-] STEP 1: Running Quality Control (FastQC)...")
    for sample in SAMPLE_INFO.index:
        r1 = f"{FASTQ_DIR}/{sample}_R1.fastq"
        r2 = f"{FASTQ_DIR}/{sample}_R2.fastq"
        
        # Output directly to QC_DIR and copy to LOGS_DIR for MultiQC processing
        cmd = f"fastqc {r1} {r2} -o {QC_DIR}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
        
    # Copy over fastqc zip outputs into logs directory for MultiQC scan path
    for file in os.listdir(QC_DIR):
        if file.endswith(".zip"):
            shutil.copy(os.path.join(QC_DIR, file), os.path.join(LOGS_DIR, file))
            
    print("[+] QC completed. Reports stored in step 1 folder.")


# =====================================================================
# STEP 2: TRIMMING & CLEANING (fastp)
# =====================================================================
def run_trimming():
    print("\n[-] STEP 2: Running Trimming and Cleaning (fastp)...")
    for sample in SAMPLE_INFO.index:
        in_r1 = f"{FASTQ_DIR}/{sample}_R1.fastq"
        in_r2 = f"{FASTQ_DIR}/{sample}_R2.fastq"
        
        out_r1 = f"{TRIM_DIR}/{sample}_trimmed_R1.fastq"
        out_r2 = f"{TRIM_DIR}/{sample}_trimmed_R2.fastq"
        html_rep = f"{TRIM_DIR}/{sample}_fastp.html"
        # Save JSON directly to the logs folder so MultiQC handles parsing
        json_rep = f"{LOGS_DIR}/{sample}_fastp.json"
        
        cmd = f"fastp -i {in_r1} -I {in_r2} -o {out_r1} -O {out_r2} -h {html_rep} -j {json_rep} --thread 4"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("[+] Trimming completed. Clean sequences ready.")


# =====================================================================
# STEP 3: ALIGNMENT / MAPPING (HISAT2)
# =====================================================================
def run_alignment():
    print("\n[-] STEP 3: Running Splice-Aware Alignment (HISAT2)...")
    index_prefix = f"{REF_DIR}/hisat2_index"
    
    if not os.path.exists(f"{index_prefix}.1.ht2"):
        print("    Building HISAT2 genomic indices...")
        build_cmd = f"hisat2-build {GENOME_FA} {index_prefix}"
        subprocess.run(build_cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
        
    for sample in SAMPLE_INFO.index:
        r1 = f"{TRIM_DIR}/{sample}_trimmed_R1.fastq"
        r2 = f"{TRIM_DIR}/{sample}_trimmed_R2.fastq"
        sam_output = f"{ALIGN_DIR}/{sample}.sam"
        log_file = f"{LOGS_DIR}/{sample}_hisat2.log"
        
        print(f"    Aligning {sample}...")
        # Pipe HISAT2 standard error (which contains alignment metrics) directly into log file
        align_cmd = f"hisat2 -x {index_prefix} -1 {r1} -2 {r2} -S {sam_output} --threads 4 2> {log_file}"
        subprocess.run(align_cmd, shell=True, check=True)
    print("[+] Mapping complete. SAM sequence alignment records generated.")


# =====================================================================
# STEP 4: QUANTIFICATION (featureCounts)
# =====================================================================
def run_quantification():
    print("\n[-] STEP 4: Running Quantification (featureCounts)...")
    sam_files = " ".join([f"{ALIGN_DIR}/{sample}.sam" for sample in SAMPLE_INFO.index])
    output_matrix_file = f"{QUANT_DIR}/raw_counts.txt"
    
    # Save stdout metrics to summary tracking path within the logs repository
    cmd = f"featureCounts -p -t exon -g gene_id -a {ANNOTATION_GTF} -o {output_matrix_file} {sam_files}"
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # featureCounts automatically drops a file matching filename.summary next to output_matrix_file
    if os.path.exists(f"{output_matrix_file}.summary"):
        shutil.copy(f"{output_matrix_file}.summary", f"{LOGS_DIR}/raw_counts.txt.summary")
        
    print("[+] Gene count aggregation matrix finished.")


# =====================================================================
# STEP 5: DE & DOWNSTREAM ANALYSIS (PyDESeq2)
# =====================================================================
def run_differential_expression():
    print("\n[-] STEP 5: Running Differential Expression Analysis (PyDESeq2)...")
    raw_file = f"{QUANT_DIR}/raw_counts.txt"
    
    df = pd.read_csv(raw_file, sep="\t", comment="#")
    df = df.set_index("GeneID")
    sample_cols = [f"{ALIGN_DIR}/{sample}.sam" for sample in SAMPLE_INFO.index]
    count_matrix = df[sample_cols]
    
    count_matrix.columns = SAMPLE_INFO.index
    count_matrix = count_matrix.T.astype(int)
    
    dds = DeseqDataSet(counts=count_matrix, metadata=SAMPLE_INFO, design_factors="condition", n_cpus=4)
    dds.deseq2()
    
    stat_res = DeseqStats(dds, contrast=["condition", "Treatment", "Control"], n_cpus=4)
    stat_res.summary()
    
    res_df = stat_res.results_df
    res_df.to_csv(f"{DE_DIR}/de_results.csv")
    
    plot_volcano(res_df)
    print(f"[+] Post-Processing Complete! Results compiled in: {DE_DIR}")


def plot_volcano(res_df):
    res_df = res_df.dropna(subset=['pvalue', 'log2FoldChange'])
    res_df['Significance'] = 'Not Significant'
    res_df.loc[(res_df['pvalue'] < 0.05) & (res_df['log2FoldChange'] > 0.5), 'Significance'] = 'Up-regulated'
    res_df.loc[(res_df['pvalue'] < 0.05) & (res_df['log2FoldChange'] < -0.5), 'Significance'] = 'Down-regulated'
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=res_df, x='log2FoldChange', y=-np.log10(res_df['pvalue'] + 1e-10), 
        hue='Significance', palette={'Not Significant': 'grey', 'Up-regulated': 'red', 'Down-regulated': 'blue'}
    )
    plt.title('Bulk RNA-Seq: Differential Gene Expression')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 P-value')
    plt.savefig(f"{DE_DIR}/volcano_plot.png", dpi=300)
    plt.close()


# =====================================================================
# STEP 6: HTML DASHBOARD AGGREGATION (MultiQC)
# =====================================================================
def generate_html_dashboard():
    print("\n[-] STEP 6: Aggregating Quality Metrics Into HTML Dashboard (MultiQC)...")
    
    # Run multiqc tool targeting everything captured inside logs sub-directory
    cmd = f"multiqc {LOGS_DIR} -o {DASHBOARD_DIR} --title 'Bulk RNA-seq Analysis Pipeline Report'"
    subprocess.run(cmd, shell=True, check=True)
    
    print(f"[+] HTML Dashboard Generated successfully! Open file://{os.path.abspath(DASHBOARD_DIR)}/multiqc_report.html in a web browser.")


# =====================================================================
# MAIN PIPELINE SYSTEM EXECUTION CONTROL
# =====================================================================
if __name__ == "__main__":
    try:
        initialize_workspace()
        run_quality_control()
        run_trimming()
        run_alignment()
        run_quantification()
        run_differential_expression()
        generate_html_dashboard()
        print("\n=======================================================")
        print("[SUCCESS] End-to-End Bulk RNA-seq Analysis Complete!")
        print(f"Check final outputs and charts under: {BASE_OUT}")
        print("=======================================================")
    except Exception as e:
        print(f"\n[CRITICAL FAILURE] Pipeline runtime exception thrown: {str(e)}")
