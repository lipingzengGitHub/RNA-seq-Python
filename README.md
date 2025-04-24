# 🧬 RNA-seq Bioinformatics Pipelines

This repository contains two Python-based RNA-seq analysis pipelines:

- **Part 1:** From raw FASTQ files to read counts using HISAT2 and HTSeq
- **Part 2:** From read counts to differentially expressed genes (DEGs) using DESeq2 via R

---

## 📌 Part 1: RNA-seq Reads Count Pipeline

**Script:** `RNA-seq-ReadsCount.py`  
This pipeline processes **paired-end RNA-seq FASTQ files** to produce read count matrices using standard bioinformatics tools.

### 🔄 Steps:

1. **Quality Control (FastQC)** – Checks raw FASTQ quality
2. **Genome Alignment (HISAT2)** – Aligns reads to the reference genome
3. **Read Counting (HTSeq-count)** – Counts reads mapped to annotated genes
4. **(Optional)** PCA clustering analysis

### 🔧 Requirements

### Conda-installed tools:
conda install -c bioconda fastqc
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda htseq

### Python packages:
pip install pandas scikit-learn matplotlib seaborn

⚠️ If in a managed system, use a virtual environment:

python3 -m venv env
source env/bin/activate
pip install pandas scikit-learn matplotlib seaborn

### Input Files
FastQ/: Folder containing paired-end .fq.gz files (e.g., sample1_1.fq.gz, sample1_2.fq.gz)
HISAT2-indexed reference genome (with .ht2 files)
Gene annotation file (GTF or GFF3)

### How to Run
python RNA-seq-ReadsCount.py -i Fastq -o hisat-count-dir -g /path/to/hisat2_index_prefix -a /path/to/annotation.gtf

### Output
FastQC_reports/: Quality reports per FASTQ
aligned_bam_files/: Aligned, sorted, and indexed BAM files
htseq_counts/: Gene-level read count files
merged_counts_matrix.csv: Combined count matrix
pca_plot.png: Optional PCA clustering plot

## 📌 Part 2: Differential Expression Analysis (DEG)
Script: run-DEG-analysis.py
This script performs differential gene expression analysis using DESeq2 via R, then generates summary files and heatmaps.

### Features:
Supports multiple test groups vs. a single control
Filters DEGs using: padj < 0.05 and |log2FC| ≥ 1
Outputs unfiltered, filtered results and heatmaps per comparison

### Input Format
python run-DEG-analysis.py counts.txt control1,control2,control3 case1_rep1,case1_rep2,case1_rep3 ...
First group is always the control. Others are test groups. Comma-separated within a group, space-separated between groups.

###  Requirements
R packages:
# From within R
install.packages("BiocManager")
BiocManager::install("DESeq2")
Or using conda:
conda install -c bioconda bioconductor-deseq2

Python environment:
conda install -c conda-forge rpy2
conda install pandas matplotlib seaborn scikit-learn
Or with pip:
pip install rpy2 pandas matplotlib seaborn scikit-learn

⚠️ Make sure Python and R versions are compatible with rpy2

### Input
counts.txt: Count matrix from Part 1 (rows = genes, columns = samples)
#Group names must match the sample headers in the count table

### How to Run
python run-DEG-analysis.py counts.txt WT_1,WT_2,WT_3 KO1_1,KO1_2,KO1_3 KO2_1,KO2_2,KO2_3

### Output
All results are saved in the DEG_results/ folder:


###File	Description
*_filtered.csv	Filtered DEGs with padj < 0.05 and `

*_unfiltered.csv	All DESeq2 results

*_heatmap.png	Heatmap of all filtered DEGs

summary_all_DEG.csv	Combined results from all group comparisons

summary_filtered_DEG.csv	Combined filtered DEGs from all comparisons

###  Notes
Supports any number of test groups, each compared independently to the control
DESeq2 handles normalization and dispersion estimation automatically
Heatmaps show expression levels of all filtered DEGs
