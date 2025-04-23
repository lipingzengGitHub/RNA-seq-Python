# Part 1: RNA-seq Reads Count Pipeline
Script: RNA-seq-ReadsCount.py

This is a Python-based pipeline for processing RNA-seq data from raw FASTQ files to read count generation using standard bioinformatics tools.

This script performs the following steps:

1. **Quality Control (FastQC)** – Performs quality check on raw FASTQ files.
2. **Genome Alignment (HISAT2)** – Aligns reads to the reference genome.
3. **Read Counting (HTSeq-count)** – Generates count data per gene using annotated features.
4. **(Optional)** PCA clustering analysis on the expression matrix.

The script is designed for **paired-end RNA-seq data**, where files end with `_1.fq.gz` and `_2.fq.gz`.

## 🔧 Requirements
Make sure the following software and Python packages are installed:

### Conda-installed tools
```bash
conda install -c bioconda fastqc
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda htseq

### Python packages
```bash
pip install pandas scikit-learn matplotlib seaborn

⚠️ You may need to use a virtual environment if you're in a restricted or externally managed system
```bash
python3 -m venv env
source env/bin/activate
pip install pandas scikit-learn matplotlib seaborn

🧬 Input Files
•	FastQ/ folder containing paired-end .fq.gz files (e.g., sample1_1.fq.gz, sample1_2.fq.gz)
•	HISAT2 indexed reference genome
•	GFF/GTF annotation file for use with HTSeq

▶️  How to Run
python RNA-seq-ReadsCount.py -i Fastq -o hisat-count-dir -g /path/to/hisat2_index_RefGenome -a /path/to/RefGenome_annotation.gff3
Arguments:
Argument	Description
-i / --input	Directory containing FASTQ files
-o / --output	Output directory for htseq-count results
-g / --genome	Path to HISAT2 genome index (without .1.ht2 etc.)
-a / --gff	GFF3 or GTF file for gene annotations

📂 Output
•	FastQC_reports/: Quality check reports for each FASTQ file
•	aligned_bam_files/: BAM files generated from HISAT2
•	hisat-count-dir/: Read count tables for each sample
•	pca_plot.png: PCA clustering of samples based on expression matrix (if used)

📣 Notes
•	Make sure the reference genome is properly indexed with hisat2-build.
•	BAM files are filtered and sorted using samtools.
•	HTSeq expects the BAM files to be sorted by position.



#
#pip install rpy2
conda install -c conda-forge rpy2

conda install -c bioconda bioconductor-deseq2




