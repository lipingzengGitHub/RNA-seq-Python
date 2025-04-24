#!/usr/bin/env python3

import os
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def print_usage_examples():
    print("""
RNA-seq Pipeline Usage Examples:
--------------------------------
To run all steps (FastQC, HISAT2 alignment, htseq-count, BAM indexing, and PCA):
  python RNA-seq-ReadsCount.py -i ./fastq -o ./output -g genome_index -a annotation.gtf --fastqc --align --count

To run only FastQC and HISAT2 alignment:
  python RNA-seq-ReadsCount.py -i ./fastq -o ./output -g genome_index --fastqc --align

To run alignment, BAM sort & index, htseq-count and PCA (skip FastQC):
  python RNA-seq-ReadsCount.py -i ./fastq -o ./output -g genome_index -a annotation.gtf --align --count

To run only htseq-count and PCA (requires existing BAM files):
  python RNA-seq-ReadsCount.py -i ./fastq -o ./output -a annotation.gtf --count
""")


def check_dependencies():
    required_tools = {
        "hisat2": "conda install -c bioconda hisat2",
        "samtools": "conda install -c bioconda samtools",
        "htseq-count": "conda install -c bioconda htseq",
        "bamCoverage": "conda install -c bioconda deeptools",
        "fastqc": "conda install -c bioconda fastqc"
    }
    for tool, install_cmd in required_tools.items():
        if not shutil.which(tool):
            print(f"Required tool '{tool}' not found.")
            print(f"Please install it with: {install_cmd}")
            exit(1)

def run_fastqc(fastq_dir, output_dir):
    print("[FastQC] Running quality control...")
    os.makedirs(output_dir, exist_ok=True)
    for f in os.listdir(fastq_dir):
        if f.endswith('.fq.gz'):
            print(f"[FastQC] Processing: {f}")
            subprocess.run(f"fastqc {os.path.join(fastq_dir, f)} -o {output_dir}", shell=True, check=True)
    print("[FastQC] Done.")

def run_hisat2(fastq_dir, output_dir, genome_index):
    print("[HISAT2] Starting alignment, sorting, indexing, and BigWig generation...")
    os.makedirs(output_dir, exist_ok=True)
    bigwig_dir = os.path.join(output_dir, "bigwig")
    os.makedirs(bigwig_dir, exist_ok=True)

    for f in os.listdir(fastq_dir):
        if f.endswith('_1.fq.gz'):
            r1 = os.path.join(fastq_dir, f)
            r2 = os.path.join(fastq_dir, f.replace('_1.fq.gz', '_2.fq.gz'))
            sample = f.replace('_1.fq.gz', '')
            sorted_bam = os.path.join(output_dir, f"{sample}_sorted.bam")
            bigwig_file = os.path.join(bigwig_dir, f"{sample}.bw")

            if os.path.exists(r2):
                print(f"[HISAT2] Aligning and sorting: {sample}")
                align_cmd = f"hisat2 -p 4 -x {genome_index} -1 {r1} -2 {r2} | samtools view -bS - | samtools sort -o {sorted_bam}"
                subprocess.run(align_cmd, shell=True, check=True)
                subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)

                print(f"[bamCoverage] Generating BigWig: {sample}.bw")
                bw_cmd = f"bamCoverage -b {sorted_bam} -o {bigwig_file} --normalizeUsing CPM"
                subprocess.run(bw_cmd, shell=True, check=True)

    print("[HISAT2] Alignment and BigWig generation completed.")

def count_reads_with_htseq(bam_dir, gff_file, output_dir):
    print("[HTSeq] Counting reads per gene...")
    os.makedirs(output_dir, exist_ok=True)
    for f in os.listdir(bam_dir):
        if f.endswith("_sorted.bam"):
            bam_path = os.path.join(bam_dir, f)
            out = os.path.join(output_dir, f"{f}_counts.txt")
            print(f"[HTSeq] Counting: {f}")
            cmd = f"htseq-count -f bam -r pos -s no -i gene_id {bam_path} {gff_file}"
            with open(out, 'w') as o:
                subprocess.run(cmd.split(), stdout=o, check=True)
    print("[HTSeq] Read counting complete.")

def merge_counts_and_plot(count_dir, output_dir):
    print("[PCA] Merging counts and generating PCA plot...")
    files = [f for f in os.listdir(count_dir) if f.endswith('_counts.txt')]
    all_counts = []

    for file in files:
        df = pd.read_csv(os.path.join(count_dir, file), sep='\t', header=None, index_col=0, comment='_')
        df = df.loc[~df.index.str.startswith('__')]
        df.columns = [file.replace('_counts.txt', '')]
        all_counts.append(df)

    merged = pd.concat(all_counts, axis=1)
    merged_file = os.path.join(output_dir, "merged_counts_matrix.csv")
    merged.to_csv(merged_file)

    data = StandardScaler().fit_transform(merged.T)
    pca = PCA(n_components=2).fit_transform(data)
    pca_df = pd.DataFrame(pca, columns=["PC1", "PC2"])
    pca_df["sample"] = merged.columns

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="sample", s=100)
    plt.title("PCA of RNA-seq Samples")
    pca_file = os.path.join(output_dir, "pca_plot.png")
    plt.savefig(pca_file)
    plt.close()
    print(f"[PCA] PCA plot saved to {pca_file}")

def main():
    parser = argparse.ArgumentParser(description="RNA-seq Pipeline: FastQC → Alignment → Counts → PCA")
    parser.add_argument("-i", "--input", required=True, help="Directory with FASTQ files")
    parser.add_argument("-o", "--output", default="RNAseq_results", help="Output directory")
    parser.add_argument("-g", "--genome", help="HISAT2 genome index prefix")
    parser.add_argument("-a", "--gff", help="Gene annotation GTF/GFF file for htseq-count")
    parser.add_argument("--fastqc", action="store_true", help="Run FastQC")
    parser.add_argument("--align", action="store_true", help="Run HISAT2 alignment")
    parser.add_argument("--count", action="store_true", help="Run htseq-count and PCA")
    args = parser.parse_args()

    check_dependencies()
    os.makedirs(args.output, exist_ok=True)

    if args.fastqc:
        run_fastqc(args.input, os.path.join(args.output, "FastQC_reports"))

    if args.align:
        if not args.genome:
            print("Error: Please provide --genome for HISAT2 alignment.")
        else:
            run_hisat2(args.input, os.path.join(args.output, "aligned_bam_files"), args.genome)

    if args.count:
        if not args.gff:
            print("Error: Please provide --gff for htseq-count.")
        else:
            count_reads_with_htseq(
                os.path.join(args.output, "aligned_bam_files"),
                args.gff,
                os.path.join(args.output, "htseq_counts")
            )
            merge_counts_and_plot(os.path.join(args.output, "htseq_counts"), args.output)

if __name__ == "__main__":
    main()

