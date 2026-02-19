import pandas as pd
import numpy as np
import sgkit as sg
import xarray as xr
import glob
import argparse

from sgkit.cohorts import _cohorts_to_array
from sgkit.stats.utils import assert_array_shape
from sgkit.utils import (
    conditional_merge_datasets,
    create_dataset,
    define_variable_if_absent,
)
from sgkit.window import has_windows, window_statistic

from sgkit import variables
from sgkit.stats.aggregation import (
    count_cohort_alleles,
    count_variant_alleles,
    individual_heterozygosity,
)

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='path to zarr dir', type=str)
parser.add_argument('-w', help='window size for df in kb', type=int)
parser.add_argument('-o', help='outpath', type=str)
args = parser.parse_args()

window_size = args.w*1000

# Functions
def read_beds(zarr_path, chrX):
    bed_path_x = zarr_path+"_X.bed"
    bed_path_all = zarr_path+".bed"
    bed_file_x = pd.read_csv(bed_path_x, sep="\t", names=["chrom", "start", "end"])
    bed_file_x = bed_file_x.loc[bed_file_x.chrom == chrX]
    bed_file_aut = pd.read_csv(bed_path_all, sep="\t", names=["chrom", "start", "end"])
    bed_file_aut = bed_file_aut.loc[bed_file_aut.chrom != chrX]
    bed_files = pd.concat([bed_file_aut, bed_file_x])
    return bed_files

def interval_creator(bed_l, window_size):
    # Input a bed file and the window size of intervals desired. Multiple chromosomes accepted.
    df_l = []
    for c in bed_l.chrom.unique():
        start_l, end_l = [], []
        b = bed_l.loc[bed_l["chrom"] == c].copy()
        b["w_s"] = b.end-b.start
        w_start = b.start.iloc[0]
        current_size = 0
        for i, j, k in zip(b.start, b.end, b.w_s):
            # Current window encapsulates the final stretch of the interval.
            if current_size + k >= window_size:
                start_l.append(w_start), end_l.append(i+(window_size-current_size))
                w_start = i+(window_size-current_size)
                # If the window still contains full intervals, contigous windows until it cant.
                for x in range((k-window_size+current_size)//window_size):
                    start_l.append(w_start), end_l.append(w_start+window_size)
                    w_start += window_size
                current_size = j-w_start
            # Current window does not encapsulate per definition, so it has to be added to current size but nothing else.
            else:
                current_size += k
        df_l.append(pd.DataFrame({"chrom": c, "interval_start": start_l, "interval_end": end_l}))
    df = pd.concat(df_l)
    df["interval_size"] = df.interval_end-df.interval_start
    return df


def non_ref_pop(ds):
    # Computes alt hom per individual in windows.
    alt_hom = ds.cohort_allele_count[:,:,1:]
    w_hom = window_statistic(
            alt_hom,
            np.sum,
            ds.window_start.values,
            ds.window_stop.values,
            dtype=np.int64,
            axis=0,
    )
    return w_hom


zarr_path = args.i

short_form = zarr_path.split("/")[-2].split("_")[0]
long_form = zarr_path.split("/")[-2]
print("Loading bed files")
bed_files = read_beds(zarr_path, long_form)
intervals_bed = interval_creator(bed_files, window_size)
# Loading the genetic data.
print("Loading genetic data")
df_l = []
ds_full = sg.load_dataset(zarr_path)
contig_IDs = pd.Series(ds_full.contig_id[:22].values).map(dict(zip(ds_full.contig_id.values, range(len(ds_full.contig_id.values))))).values

df_l = []
for c, c_ID in zip(ds_full.contig_id[:1].values, contig_IDs):
    print(c)
    ds = ds_full.sel(variants=((ds_full.variant_contig == c_ID).compute()))
    c_intervals = intervals_bed.loc[intervals_bed.chrom == c]
    if len(ds.variants) < 100:
        print("Skipping due to too few variants")
        continue
    var_chunk = min(50, 1+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = ds["samples"]
    # Subsetting and windowing the sgkit dataset.
    # The rechunking handles what otherwise would cause an error.
    ds["call_genotype"] = ds["call_genotype"].clip(0)
    ds = ds.sel(contigs=[ds.variant_contig[0].values])
    ds["interval_contig_name"] = (["intervals"], c_intervals.chrom)
    ds["interval_start"] = (["intervals"], c_intervals.interval_start)
    ds["interval_stop"] = (["intervals"], c_intervals.interval_end)
    ds = sg.window_by_interval(ds)
    ds = (sg.diversity(ds.chunk({"variants": len(ds.variants)//var_chunk}))) # Create at most 50 chunks
    ds["non_ref_count"] = (("windows", "cohorts"), non_ref_pop(ds))
    for i in range(len(ds.sample_id)):
        df_sub = pd.DataFrame({"het": ds.stat_diversity[:,i], "alt_hom": ds.alt_hom[:,i],
        "variant_count": ds.window_stop-ds.window_start, "GVCF_ID": ds.sample_id[i].values})
        df_sub["window_start"] = list(range(0, len(ds.window_start)*window_size, window_size))
        df_sub["chrom"] = c
        df_l.append(df_sub)

ds_X = sg.load_dataset(zarr_path+"_X")
# Human chrX
contig_IDs = pd.Series(ds_full.contig_id[22:23].values).map(dict(zip(ds_full.contig_id.values, range(len(ds_full.contig_id.values))))).values

for c, c_ID in zip(ds_full.contig_id[22:23].values, contig_IDs):
    print(c)
    ds = ds_X.sel(variants=((ds_X.variant_contig == c_ID).compute()))
    c_intervals = intervals_bed.loc[intervals_bed.chrom == c]
    if len(ds.variants) < 100:
        print("Skipping due to too few variants")
        continue
    var_chunk = min(50, 1+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = ds["samples"]
    # Subsetting and windowing the sgkit dataset.
    # The rechunking handles what otherwise would cause an error.
    ds["call_genotype"] = ds["call_genotype"].clip(0)
    ds = ds.sel(contigs=[ds.variant_contig[0].values])
    ds["interval_contig_name"] = (["intervals"], c_intervals.chrom)
    ds["interval_start"] = (["intervals"], c_intervals.interval_start)
    ds["interval_stop"] = (["intervals"], c_intervals.interval_end)
    ds = sg.window_by_interval(ds)
    ds = (sg.diversity(ds.chunk({"variants": len(ds.variants)//var_chunk}))) # Create at most 50 chunks
    ds["non_ref_count"] = (("windows", "cohorts"), non_ref_pop(ds))
    for i in range(len(ds.sample_id)):
        df_sub = pd.DataFrame({"het": ds.stat_diversity[:,i], "alt_hom": ds.alt_hom[:,i],
        "variant_count": ds.window_stop-ds.window_start, "GVCF_ID": ds.sample_id[i].values})
        df_sub["window_start"] = list(range(0, len(ds.window_start)*window_size, window_size))
        df_sub["chrom"] = c
        df_l.append(df_sub)

df_het = pd.concat(df_l)
intervals_callable = pos_windows(bed_files, window_size, bed_files["chrom"].unique())
output_df = pd.merge(df_het, intervals_callable, on=["chrom", "window_start"])
output_df["chr_type"] = ["chrX" if x == 'NC_060947.1' else "aut" for x in output_df["chrom"]]
output_df["species"] = long_form
output_df.to_csv("{}{}_{}kb_het_hom.txt".format(args.o, long_form, args.w), sep="\t", index=False)
