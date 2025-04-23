def print_usage_examples():
    print("""
Pipeline Usage Examples:
--------------------------------

python run-DEG-analysis.py counts.txt control_1,control_2,control_3 treatmentA_1,treatmentA_2,treatmentA_3 treatmentB_1,treatmentB_2,treatmentB_3


""")

import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri, r, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

# Activate automatic conversion between R and pandas
pandas2ri.activate()
base = importr('base')
utils = importr('utils')
deseq2 = importr('DESeq2')

# Create result directory
os.makedirs("DEG_results", exist_ok=True)

def run_deseq(counts, col_data, control_name, case_name, out_prefix):
    dds = deseq2.DESeqDataSetFromMatrix(countData=counts,
                                         colData=col_data,
                                         design=Formula('~ group'))
    dds = deseq2.DESeq(dds)

    res = deseq2.results(dds, contrast=robjects.StrVector(['group', case_name, control_name]))
    res_df_r = base.as_data_frame(res)
    res_df = pandas2ri.rpy2py(res_df_r)

    res_df.index.name = 'gene'
    res_df = res_df.sort_values('padj')

    # Save unfiltered results
    res_df.to_csv(f"DEG_results/{out_prefix}_unfiltered.csv")

    # Filter DEGs: padj < 0.05 and |log2FC| ≥ 1
    res_filtered = res_df[(res_df['padj'] < 0.05) & (abs(res_df['log2FoldChange']) >= 1)]
    res_filtered.to_csv(f"DEG_results/{out_prefix}_filtered.csv")

    return res_df, res_filtered

def generate_heatmap(counts_df, genes, out_file):
    if len(genes) == 0:
        print(f"[Heatmap] No DEGs to plot for {out_file}")
        return
    data = counts_df.loc[genes]
    sns.clustermap(data, standard_scale=0, cmap="vlag", figsize=(10, 8))
    plt.savefig(out_file)
    plt.close()
    print(f"[Heatmap] Saved to {out_file}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python run-DEG-analysis.py counts.txt control_rep1,control_rep2,... case1_rep1,... case2_rep1,...")
        sys.exit(1)

    count_file = sys.argv[1]
    group_inputs = sys.argv[2:]

    # Sample group parsing
    groups = [grp.split(',') for grp in group_inputs]
    group_labels = [f"group{i+1}" for i in range(len(groups))]

    # Read count table
    count_df = pd.read_csv(count_file, sep='\t', index_col=0)

    summary_all = []
    summary_filtered = []

    control_name = group_labels[0]
    control_samples = groups[0]

    for i in range(1, len(groups)):
        case_name = group_labels[i]
        case_samples = groups[i]

        samples = control_samples + case_samples
        count_subset = count_df[samples]

        group_values = [control_name]*len(control_samples) + [case_name]*len(case_samples)
        col_data = pd.DataFrame({'group': pd.Categorical(group_values)}, index=samples)

        r_count = pandas2ri.py2rpy(count_subset.astype(int))
        r_coldata = pandas2ri.py2rpy(col_data)

        out_prefix = f"{case_name}_vs_{control_name}"
        print(f"[DESeq2] Running: {out_prefix}")
        res_all, res_filtered = run_deseq(r_count, r_coldata, control_name, case_name, out_prefix)

        summary_all.append(res_all.rename(columns={
            "log2FoldChange": f"log2FC_{out_prefix}",
            "padj": f"padj_{out_prefix}"
        }))
        summary_filtered.append(res_filtered)

        # Plot heatmap of all filtered DEGs (not top 50)
        generate_heatmap(count_subset, res_filtered.index, f"DEG_results/{out_prefix}_heatmap.png")

    # Save summary files
    if summary_all:
        merged_all = pd.concat(summary_all, axis=1)
        merged_all.to_csv("DEG_results/summary_all_DEG.csv")
    if summary_filtered:
        all_filtered_df = pd.concat(summary_filtered).drop_duplicates()
        all_filtered_df.to_csv("DEG_results/summary_filtered_DEG.csv")

if __name__ == "__main__":
    main()

