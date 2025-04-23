import os
import gzip
import random

def generate_dummy_fastq(output_dir, sample_prefix, num_reads=1000):
    """
    Generate paired-end dummy FASTQ files for testing.
    Each read will be a random sequence of 50 bases with corresponding quality scores.
    """
    os.makedirs(output_dir, exist_ok=True)
    bases = ['A', 'T', 'C', 'G']
    
    for pair in ['_1', '_2']:
        file_path = os.path.join(output_dir, f"{sample_prefix}{pair}.fq.gz")
        with gzip.open(file_path, 'wt') as f:
            for i in range(num_reads):
                seq = ''.join(random.choices(bases, k=50))
                qual = ''.join(chr(random.randint(33, 73)) for _ in range(50))
                f.write(f"@{sample_prefix}_read{i}{pair}\n{seq}\n+\n{qual}\n")

def generate_test_fastq_set():
    """
    Generate a test dataset with 3 paired-end samples.
    """
    test_dir = "test_fastq"
    samples = ["sample1", "sample2", "sample3"]

    for sample in samples:
        generate_dummy_fastq(test_dir, sample)
    
    print(f"Dummy FASTQ files generated in ./{test_dir}")

generate_test_fastq_set()


