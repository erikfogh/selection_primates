import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', help='path to file', type=str)
parser.add_argument('-g', help='path to feature bed', type=str)
parser.add_argument('-w', help='interval grouping size', type=int)
args = parser.parse_args()

interval_size = args.w*1000
pi_df = pd.read_csv(args.i, sep="\t")
CAT_Genes = pd.read_csv(args.g, sep="\t")
protein_coding = CAT_Genes.loc[CAT_Genes.geneType == "protein_coding"]
chr_names = ['NC_060925.1', 'NC_060926.1', 'NC_060927.1', 'NC_060928.1',
       'NC_060929.1', 'NC_060930.1', 'NC_060931.1', 'NC_060932.1',
       'NC_060933.1', 'NC_060934.1', 'NC_060935.1', 'NC_060936.1',
       'NC_060937.1', 'NC_060938.1', 'NC_060939.1', 'NC_060940.1',
       'NC_060941.1', 'NC_060942.1', 'NC_060943.1', 'NC_060944.1',
       'NC_060945.1', 'NC_060946.1', 'NC_060947.1']
chr_numbers = ["chrX" if i == 22 else "chr{}".format(i+1) for i in range(len(chr_names))]

df_l = []
for c_n, c_i in zip(chr_names[:], chr_numbers[:]):
    print(c_i)
    chr_genes = protein_coding.loc[protein_coding["#chrom"] == c_i].copy()
    chr_df = pi_df.loc[pi_df.chrom == c_n].copy()
    if len(chr_df) == 0:
        print("Skipping due to no data")
        continue
    chr_df["window_size"] = chr_df.window_end-chr_df.window_start
    chr_df["midpoint"] = (chr_df.window_end-chr_df.window_start)/2+chr_df.window_start
    # This windowing procedure only makes sense if windows larger than 100kb are filtered out, but I will do that.
    chr_df["window_{}kb".format(args.w)] = (pd.cut(chr_df["midpoint"],
                                     np.arange(0, chr_df.window_end.max()+interval_size, interval_size),
                                     labels=False).values)*interval_size
    genes_l, max_genes_l, largest_gene_l, largest_gene_cov_l = [], [], [], []
    for w_s, w_e in zip(chr_df.window_start, chr_df.window_end):
        overlapping_bed = chr_genes.loc[(chr_genes.thickStart <= w_e) & (chr_genes.thickEnd >= w_s)].copy()
        genes = overlapping_bed.name2.unique()
        if len(genes) >= 1:
            overlapping_bed["clipped_start"] = overlapping_bed.thickStart.clip(w_s, w_e)
            overlapping_bed["clipped_end"] = overlapping_bed.thickEnd.clip(w_s, w_e)
            overlapping_bed["feature_length"] = overlapping_bed.clipped_end-overlapping_bed.clipped_start
            sorted_max = overlapping_bed.groupby(["name2"])[["feature_length"]].max().sort_values(by="feature_length",
                                                                                     ascending=False).reset_index()
            genes_l.append(sorted_max.name2.values)
            max_genes_l.append(len(sorted_max.loc[sorted_max.feature_length >= w_e-w_s]))
            largest_gene_l.append(sorted_max.iloc[0].name2)
            largest_gene_cov_l.append(sorted_max.iloc[0].feature_length/(w_e-w_s))
        else:
            genes_l.append(np.nan)
            max_genes_l.append(0)
            largest_gene_l.append(np.nan)
            largest_gene_cov_l.append(0)
    chr_df["genes"] = genes_l
    chr_df["max_genes"] = max_genes_l
    chr_df["largest_gene"] = largest_gene_l
    chr_df["largest_gene_cov"] = largest_gene_cov_l
    df_l.append(chr_df)
annotated_df = pd.concat(df_l)
annotated_df.to_csv(args.i[:-4]+"_annotated.txt", sep="\t", index=False)
