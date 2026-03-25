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


def fix_intervals(ds_cohort_div, c_intervals, window_contig):
    k = 0
    new_start, new_end, contig_l = [], [], []
    variant_pos_vals = ds_cohort_div.variant_position[ds_cohort_div.window_start.values].values
    for i, j, l in zip(ds_cohort_div.window_start.values, ds_cohort_div.window_stop.values, variant_pos_vals):
        while c_intervals.iloc[k].interval_end < l:
            #print(c_intervals.iloc[k].interval_end, ds_cohort_div.variant_position[i].values, "fail")
            new_start.append(i), new_end.append(i), contig_l.append(window_contig)
            k += 1
        #print(c_intervals.iloc[k].interval_end, ds_cohort_div.variant_position[i].values, "pass")
        new_start.append(i), new_end.append(j), contig_l.append(window_contig)
        k += 1
    ds_cohort_div = ds_cohort_div.drop_dims("windows")
    ds_cohort_div["window_contig"] = ("windows", contig_l)
    ds_cohort_div["window_start"] = ("windows", new_start)
    ds_cohort_div["window_stop"] = ("windows", new_end)
    ds_cohort_div = ds_cohort_div.assign_attrs(units="windows")
    return ds_cohort_div


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


parser = argparse.ArgumentParser()
parser.add_argument('-i', help='path to zarr dir', type=str)
parser.add_argument('-w', help='window size for df in kb', type=int)
parser.add_argument('-o', help='outpath', type=str)
args = parser.parse_args()

window_size = args.w*1000
zarr_path = args.i
short_form = zarr_path.split("/")[-2].split("_")[0]
long_form = zarr_path.split("/")[-2]
print("Loading bed files")
bed_files = read_beds(zarr_path, 'NC_060947.1')
intervals_bed = interval_creator(bed_files, window_size)
# Loading the genetic data.
print("Loading genetic data")
ds_full = sg.load_dataset(zarr_path)
df_l = []

contig_IDs = pd.Series(ds_full.contig_id[:22].values).map(dict(zip(ds_full.contig_id.values, range(len(ds_full.contig_id.values))))).values
for c, c_ID in zip(ds_full.contig_id[:22].values, contig_IDs):
    print(c)
    ds = ds_full.sel(variants=((ds_full.variant_contig == c_ID).compute()))
    c_intervals = intervals_bed.loc[intervals_bed.chrom == c]
    if len(ds.variants) < 1000:
        print("Skipping due to too few variants")
        continue
    var_chunk = min(50, 1+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = [0]*len(ds["samples"])
    # Subsetting and windowing the sgkit dataset.
    # The rechunking handles what otherwise would cause an error.
    ds = ds.sel(contigs=[ds.variant_contig[0].values])
    ds["interval_contig_name"] = (["intervals"], c_intervals.chrom)
    ds["interval_start"] = (["intervals"], c_intervals.interval_start)
    ds["interval_stop"] = (["intervals"], c_intervals.interval_end)
    ds_cohort = sg.count_cohort_alleles(ds)
    ds_cohort_div = sg.window_by_interval(ds_cohort)
    #Intervals with no variant sites are not kept with the interval method, this adds them
    ds_cohort_div = fix_intervals(ds_cohort_div, c_intervals, c_ID)
    ds_cohort_div["alt_hom"] = (("windows", "cohorts"), non_ref_pop(ds_cohort_div)[:,:,0])
    chrom_df = pd.DataFrame({"window_start": c_intervals.interval_start,
                            "window_end": c_intervals.interval_end})
    chrom_df["divergence"] = ds_cohort_div.alt_hom/len(ds["samples"])
    ds_cohort_pi = ds_cohort.sel(variants = ((ds_cohort.cohort_allele_count[:,0] >= 1).sum(axis=1) >= 2).compute())
    if len(ds_cohort_pi.variants) < 10:
        chrom_df["pi"] = np.nan
    else:
        ds_cohort_pi = sg.window_by_interval(ds_cohort_pi) 
        ds_cohort_pi = fix_intervals(ds_cohort_pi, c_intervals, c_ID)
        ds_cohort_pi = (sg.diversity(ds_cohort_pi))
        chrom_df["pi"] = pd.Series(ds_cohort_pi.stat_diversity[:,0])
    chrom_df["chrom"] = c
    df_l.append(chrom_df)

ds_X = sg.load_dataset(zarr_path+"_X")
# Human chrX
contig_IDs = pd.Series(ds_full.contig_id[22:23].values).map(dict(zip(ds_full.contig_id.values, range(len(ds_full.contig_id.values))))).values

for c, c_ID in zip(ds_full.contig_id[22:23].values, contig_IDs):
    print(c)
    ds = ds_X.sel(variants=((ds_X.variant_contig == c_ID).compute()))
    c_intervals = intervals_bed.loc[intervals_bed.chrom == c]
    if len(ds.variants) < 1000:
        print("Skipping due to too few variants")
        continue
    var_chunk = min(50, 1+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = [0]*len(ds["samples"])
    # Subsetting and windowing the sgkit dataset.
    # The rechunking handles what otherwise would cause an error.
    ds = ds.sel(contigs=[ds.variant_contig[0].values])
    ds["interval_contig_name"] = (["intervals"], c_intervals.chrom)
    ds["interval_start"] = (["intervals"], c_intervals.interval_start)
    ds["interval_stop"] = (["intervals"], c_intervals.interval_end)
    ds_cohort = sg.count_cohort_alleles(ds)
    ds_cohort_div = sg.window_by_interval(ds_cohort)
    ds_cohort_div = fix_intervals(ds_cohort_div, c_intervals, c_ID)
    ds_cohort_div["alt_hom"] = (("windows", "cohorts"), non_ref_pop(ds_cohort_div)[:,:,0])
    chrom_df = pd.DataFrame({"window_start": c_intervals.interval_start,
                            "window_end": c_intervals.interval_end})
    chrom_df["divergence"] = ds_cohort_div.alt_hom/len(ds["samples"])
    ds_cohort_pi = ds_cohort.sel(variants = ((ds_cohort.cohort_allele_count[:,0] >= 1).sum(axis=1) >= 2).compute())
    if len(ds_cohort_pi.variants) < 10:
        chrom_df["pi"] = np.nan
    else:
        ds_cohort_pi = sg.window_by_interval(ds_cohort_pi) 
        ds_cohort_pi = fix_intervals(ds_cohort_pi, c_intervals, c_ID)
        ds_cohort_pi = (sg.diversity(ds_cohort_pi))
        chrom_df["pi"] = pd.Series(ds_cohort_pi.stat_diversity[:,0])
    chrom_df["chrom"] = c
    df_l.append(chrom_df)

output_df = pd.concat(df_l)
output_df["chr_type"] = ["chrX" if x == 'NC_060947.1' else "aut" for x in output_df["chrom"]]
output_df["species"] = long_form
output_df.to_csv("{}{}_{}kb_pi.txt".format(args.o, long_form, args.w), sep="\t", index=False)
