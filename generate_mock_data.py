import os

# Create directory structure
os.makedirs("./fastq", exist_ok=True)
os.makedirs("./reference", exist_ok=True)

# 1. Create a tiny mock reference genome (FASTA)
with open("./reference/genome.fa", "w") as f:
    f.write(">chr1\n" + "ATGC" * 25000 + "\n")  # A simple 10kb genome string

# 2. Create mock gene annotations (GTF file format required by featureCounts)
with open("./reference/annotations.gtf", "w") as f:
    f.write('chr1\tGenbank\tgene\t100\t2000\t.\t+\t.\tgene_id "gene_A"; transcript_id "tx_A";\n')
    f.write('chr1\tGenbank\texon\t100\t2000\t.\t+\t.\tgene_id "gene_A"; transcript_id "tx_A";\n')
    f.write('chr1\tGenbank\tgene\t3000\t5000\t.\t+\t.\tgene_id "gene_B"; transcript_id "tx_B";\n')
    f.write('chr1\tGenbank\texon\t3000\t5000\t.\t+\t.\tgene_id "gene_B"; transcript_id "tx_B";\n')

# 3. Create mock paired-end sequencing reads (FASTQ format)
samples = ["sample1", "sample2", "sample3", "sample4"]
for sample in samples:
    for read in ["R1", "R2"]:
        with open(f"./fastq/{sample}_{read}.fastq", "w") as f:
            for i in range(100):  # Generate 100 mock reads per file
                f.write(f"@{sample}_read_{i}\n")
                f.write("ATGC" * 25 + "\n")  # 100bp sequence read
                f.write("+\n")
                f.write("I" * 100 + "\n")    # High quality scores (Phred+33)

print("[+] Mock dataset generated successfully in ./fastq and ./reference!")
