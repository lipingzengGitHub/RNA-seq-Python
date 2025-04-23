import os
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

def run_fastqc(fastq_dir, output_dir):
    print("[FastQC] Starting quality control...")
    os.makedirs(output_dir, exist_ok=True)
    for f in os.listdir(fastq_dir):
        if f.endswith('.fq.gz'):
            print(f"[FastQC] Processing: {f}")
            subprocess.run(f"fastqc {os.path.join(fastq_dir, f)} -o {output_dir}", shell=True, check=True)

def run_hisat2(fastq_dir, output_dir, genome_index):
    print("[HISAT2] Starting alignment + sorting + indexing + BigWig...")
    os.makedirs(output_dir, exist_ok=True)
    for f in os.listdir(fastq_dir):
        if f.endswith('_1.fq.gz'):
            r1 = os.path.join(fastq_dir, f)
            r2 = os.path.join(fastq_dir, f.replace('_1.fq.gz', '_2.fq.gz'))
            sample = f.replace('_1.fq.gz', '')
            sorted_bam = os.path.join(output_dir, f"{sample}_sorted.bam")
            if os.path.exists(r2):
                print(f"[HISAT2] Aligning and sorting: {sample}")
                cmd = f"hisat2 -p 4 -x {genome_index} -1 {r1} -2 {r2} | samtools view -bS - | samtools sort -o {sorted_bam}"
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
                os.makedirs(os.path.join(output_dir, "bigwig"), exist_ok=True)
                subprocess.run(f"bamCoverage -b {sorted_bam} -o {output_dir}/bigwig/{sample}.bw --normalizeUsing CPM", shell=True, check=True)
                print(f"[HISAT2] Done: {sample}")
            else:
                print(f"[HISAT2] Missing R2 for {f}, skipped.")

def count_reads_with_htseq(bam_dir, gff_file, output_dir):
    print("[HTSeq] Counting reads...")
    os.makedirs(output_dir, exist_ok=True)
    for f in os.listdir(bam_dir):
        if f.endswith(".bam"):
            bam_path = os.path.join(bam_dir, f)
            output_path = os.path.join(output_dir, f"{f}_counts.txt")
            print(f"[HTSeq] Processing: {f}")
            cmd = f"htseq-count -f bam -r pos -s no -i gene_id {bam_path} {gff_file}"
            with open(output_path, 'w') as out:
                subprocess.run(cmd.split(), stdout=out, check=True)

def merge_counts_and_plot(count_dir, output_dir):
    print("[Merge] Merging count tables and PCA...")
    files = [f for f in os.listdir(count_dir) if f.endswith("_counts.txt")]
    all_counts = []
    for file in files:
        df = pd.read_csv(os.path.join(count_dir, file), sep='\t', header=None, index_col=0, comment='_')
        df = df[~df.index.str.startswith("__")]
        df.columns = [file.replace("_sorted.bam_counts.txt", "").replace("_counts.txt", "")]
        all_counts.append(df)
    merged = pd.concat(all_counts, axis=1)
    merged.to_csv(os.path.join(output_dir, "merged_counts_matrix.csv"))
    # PCA
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(merged.T)
    pca = PCA(n_components=2).fit_transform(data_scaled)
    df_pca = pd.DataFrame(pca, columns=["PC1", "PC2"])
    df_pca["sample"] = merged.columns
    plt.figure(figsize=(8,6))
    sns.scatterplot(x="PC1", y="PC2", data=df_pca, hue="sample")
    plt.title("PCA of RNA-seq Samples")
    plt.savefig(os.path.join(output_dir, "pca_plot.png"))
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input FASTQ directory")
    parser.add_argument("-o", "--output", default="RNAseq_results", help="Output directory")
    parser.add_argument("-g", "--genome", help="HISAT2 index prefix")
    parser.add_argument("-a", "--gff", help="Annotation GFF/GTF file")
    parser.add_argument("--fastqc", action="store_true", help="Run FastQC")
    parser.add_argument("--align", action="store_true", help="Run HISAT2 + sort + index + BigWig")
    parser.add_argument("--count", action="store_true", help="Run htseq-count and PCA")
    args = parser.parse_args()

    print_usage_examples()

    os.makedirs(args.output, exist_ok=True)
    if args.fastqc:
        run_fastqc(args.input, os.path.join(args.output, "FastQC_reports"))
    if args.align:
        if not args.genome:
            print("[ERROR] Genome index (-g) is required for alignment.")
        else:
            run_hisat2(args.input, os.path.join(args.output, "aligned_bam_files"), args.genome)
    if args.count:
        if not args.gff:
            print("[ERROR] Annotation file (-a) is required for htseq-count.")
        else:
            bam_dir = os.path.join(args.output, "aligned_bam_files")
            count_dir = os.path.join(args.output, "htseq_counts")
            count_reads_with_htseq(bam_dir, args.gff, count_dir)
            merge_counts_and_plot(count_dir, args.output)

if __name__ == "__main__":
    main()


