import pandas as pd
import numpy as np
import os

# Create directory for test data
os.makedirs("test_data_DEG", exist_ok=True)

# Gene list and sample names
genes = [f"Gene{i}" for i in range(1, 101)]
control_samples = ["1-1-1", "1-1-2", "1-1-3"]
case1_samples = ["1-2-1", "1-2-2", "1-2-3"]
case2_samples = ["1-3-1", "1-3-2", "1-3-3"]
samples = control_samples + case1_samples + case2_samples

# Simulate expression counts: baseline for control, elevated in cases
np.random.seed(42)
counts = pd.DataFrame(index=genes, columns=samples)
for gene in genes:
    base = np.random.poisson(lam=100)
    counts.loc[gene, control_samples] = np.random.poisson(lam=base, size=3)
    counts.loc[gene, case1_samples] = np.random.poisson(lam=base*2, size=3)
    counts.loc[gene, case2_samples] = np.random.poisson(lam=base*0.5, size=3)

# Save to counts.txt in tab-separated format
counts = counts.astype(int)
counts.to_csv("test_data_DEG/counts.txt", sep='\t')

