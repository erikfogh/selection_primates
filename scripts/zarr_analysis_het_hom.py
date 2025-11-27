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
parser.add_argument('-m', help='path to metadata', type=str)
parser.add_argument('-w', help='window size for df in kb', type=int)
parser.add_argument('-o', help='outpath', type=str)
args = parser.parse_args()

window_size = args.w*1000

# Functions

def read_beds(long_form):
    bed_path_x = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/{}/filteredVCF/pos_bed_cov_based/{}_batch*_fploidy2_mploidy1.bed".format(long_form, long_form)
    bed_path_all = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/{}/filteredVCF/pos_bed_cov_based/{}_batch*_fploidy2_mploidy2.bed".format(long_form, long_form)
    bed_l = []
    for b in glob.glob(bed_path_all):
        bed_file = pd.read_csv(b, sep="\t", names=["chrom", "start", "end"])
        bed_l.append(bed_file)
    bed_files = pd.concat(bed_l)
    bed_l = []
    for b in glob.glob(bed_path_x):
        #print(b)
        bed_file = pd.read_csv(b, sep="\t", names=["chrom", "start", "end"])
        bed_l.append(bed_file)
    if len(bed_l) > 0:
        bed_x = pd.concat(bed_l)
        bed_files = bed_files.loc[~(bed_files.chrom.isin(bed_x.chrom.unique()))]
        bed_files = pd.concat([bed_files, bed_x]).sort_values(by=["chrom", "start", "end"])
    return bed_files


def pos_windows(bed_l, window_size, chrom_order):
    # Input a bed file and the window size of intervals desired. Multiple chromosomes accepted.
    # It has to be sorted.
    df_l = []
    for c in chrom_order:
        #print(c)
        frac_l = []
        b = bed_l.loc[bed_l["chrom"] == c].copy()
        b["w_s"] = b.end-b.start
        w_start = b.start.iloc[0]
        current_pos, callable_bases = 0, 0
        for i, j, k in zip(b.start, b.end, b.w_s):
            # Nothing called in the current window under investigation.
            while i-window_size >= current_pos:
                frac_l.append(callable_bases/window_size)
                callable_bases = 0
                current_pos += window_size
            # Window starts in current. We know this is true because of the previous while loop.
            callable_bases += min(k, current_pos+window_size-i)
            # Everything called in current.
            while j-window_size >= current_pos:
                frac_l.append(callable_bases/window_size)
                callable_bases = 0
                current_pos += window_size
                if j-window_size >= current_pos:
                    callable_bases += window_size
                else:
                # Window stops in current. Again, know this is true.
                    callable_bases += j-current_pos
        # Last window.
        frac_l.append(callable_bases/(window_size))
        df_l.append(pd.DataFrame({"chrom": c, "window_start": list(range(0, len(frac_l)*window_size, window_size)),
                                  "window_end": list(range(window_size, (len(frac_l)+1)*window_size, window_size)),
                                  "callable_frac": frac_l}))
    return pd.concat(df_l)


def alt_hom(ds):
    # Computes alt hom per individual in windows.
    alt_hom = ((ds.call_genotype[:,:,0] == ds.call_genotype[:,:,-1]) & (ds.call_genotype[:,:,0] >= 1))
    w_hom = window_statistic(
            alt_hom,
            np.sum,
            ds.window_start.values,
            ds.window_stop.values,
            dtype=np.int64,
            axis=0,
    )
    return w_hom


print("Loading metadata")
zarr_path = args.i
short_form = zarr_path.split("/")[-1].split("_")[0]
long_form = zarr_path.split("/")[-1]
# Loading the various metadata files. Metadata, contig information, callability bed.
metadata_path = args.m
metadata_df = pd.read_csv(metadata_path+"{}_individuals.txt".format(short_form), sep="\t")
metadata_df["SEX_I"] = [0 if x == "F" else 1 for x in metadata_df.GENETIC_SEX]
regions_df = pd.read_csv(metadata_path+"{}_regions_and_batches.txt".format(short_form), sep="\t")
regions_df["LENGTH"] = regions_df["END"]-regions_df["START"]
regions_df["chr_type"] = ["chrX" if x == 2 and y == 1 else "aut" for x, y in zip(regions_df.FEMALE_PLOIDY, regions_df.MALE_PLOIDY)]
large_contigs = regions_df.loc[(regions_df.LENGTH >= 1000000) & (regions_df.FEMALE_PLOIDY == 2)].CONTIG_ID.unique()
large_x = regions_df.loc[(regions_df.LENGTH >= 1000000) & (regions_df.FEMALE_PLOIDY == 2) &
                        (regions_df.MALE_PLOIDY == 1)].CONTIG_ID
print("Loading bed files")
bed_files = read_beds(long_form)
# Loading the genetic data.
print("Loading genetic data")
df_l = []
for c in glob.glob(zarr_path+"/*/"):
    print(c)
    ds = sg.load_dataset(c[:-1])
    if len(ds.variants) < 10:
        print("Skipping due to too few variants")
    var_chunk = min(50, 1+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = ds["samples"]
    # Subsetting and windowing the sgkit dataset.
    # The rechunking handles what otherwise would cause an error.
    ds["call_genotype"] = ds["call_genotype"].clip(0)
    ds = ds.sel(contigs=[ds.variant_contig[0].values])
    ds = sg.window_by_position(ds, size=window_size)
    ds = (sg.diversity(ds.chunk({"variants": len(ds.variants)//var_chunk}))) # Create at most 50 chunks
    ds["alt_hom"] = (("windows", "cohorts"), alt_hom(ds))
    for i in range(len(ds.sample_id)):
        df_sub = pd.DataFrame({"het": ds.stat_diversity[:,i], "alt_hom": ds.alt_hom[:,i],
        "variant_count": ds.window_stop-ds.window_start, "GVCF_ID": ds.sample_id[i].values})
        df_sub["window_start"] = list(range(0, len(ds.window_start)*window_size, window_size))
        df_sub["chrom"] = c.split("/")[-1]
        df_l.append(df_sub)
df_het = pd.concat(df_l)
bed_files = read_beds(long_form)
intervals_callable = pos_windows(bed_files, window_size, bed_files["chrom"].unique())
output_df = pd.merge(df_het, intervals_callable, on=["chrom", "window_start"])
output_df["chr_type"] = output_df["chrom"].map(dict(zip(regions_df.CONTIG_ID, regions_df.chr_type)))
output_df["species"] = long_form
output_df.to_csv("{}{}_{}kb_het_hom.txt".format(args.o, long_form, args.w), sep="\t", index=False)
