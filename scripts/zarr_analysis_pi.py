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

from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering


parser = argparse.ArgumentParser()
parser.add_argument('-i', help='path to zarr dir', type=str)
parser.add_argument('-m', help='path to metadata', type=str)
parser.add_argument('-w', help='window size for df in kb', type=int)
parser.add_argument('-o', help='outpath', type=str)
args = parser.parse_args()

window_size = args.w*1000
bed_base = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/"

# Functions

def read_beds(long_form):
    bed_path_x = bed_base+"{}/filteredVCF/all_samples/pos_bed_cov_based/{}_batch*_fploidy2_mploidy1.bed".format(long_form, long_form)
    bed_path_all = bed_base+"{}/filteredVCF/all_samples/pos_bed_cov_based/{}_batch*_fploidy2_mploidy2.bed".format(long_form, long_form)
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


print("Loading metadata")
zarr_path = args.i
short_form = zarr_path.split("/")[-2].split("_")[0]
long_form = zarr_path.split("/")[-2]
# Loading the various metadata files. Metadata, contig information, callability bed.
metadata_path = args.m
metadata_df = pd.read_csv(metadata_path+"{}_individuals.txt".format(short_form), sep="\t")
metadata_df["SEX_I"] = [0 if x == "F" else 1 for x in metadata_df.GENETIC_SEX]
regions_df = pd.read_csv(metadata_path+"{}_regions_and_batches.txt".format(short_form), sep="\t")
regions_df["LENGTH"] = regions_df["END"]-regions_df["START"]
regions_df["chr_type"] = ["chrX" if x == 2 and y == 1 else "aut" for x, y in zip(regions_df.FEMALE_PLOIDY, regions_df.MALE_PLOIDY)]
large_contigs = regions_df.loc[(regions_df.LENGTH >= 1000000) & (regions_df.FEMALE_PLOIDY == 2)].CONTIG_ID.unique()
large_aut = regions_df.loc[(regions_df.LENGTH >= 1000000) & (regions_df.FEMALE_PLOIDY == 2) &
                        (regions_df.MALE_PLOIDY == 2)].CONTIG_ID
large_x = regions_df.loc[(regions_df.LENGTH >= 1000000) & (regions_df.FEMALE_PLOIDY == 2) &
                        (regions_df.MALE_PLOIDY == 1)].CONTIG_ID
print("Loading bed files")
bed_files = read_beds(long_form)
# Loading the genetic data.
print("Loading genetic data")
ds_full = sg.load_dataset(zarr_path)

kept_contigs =  [x for x in ds_full.contig_id.values if (x == large_contigs).any()]
contig_IDs = pd.Series(kept_contigs).map(dict(zip(ds_full.contig_id.values, range(len(ds_full.contig_id.values))))).values


df_l = []
for c, c_ID in zip(kept_contigs, contig_IDs):
    print(c)
    ds = ds_full.sel(variants=((ds_full.variant_contig == c_ID).compute()), contigs=[c_ID])
    if len(ds.variants) < 100:
        print("Skipping due to too few variants")
        continue
    # Some chunking works as too many chunks can cause problems.
    var_chunk = min(50, 5+len(ds.variants)//1000000)
    print(len(ds.variants), var_chunk)
    ds["sample_cohort"] = [0]*len(ds["samples"])
    ds_cohort = sg.count_cohort_alleles(ds)
    # Counting non ref positions
    ds_cohort_div = sg.window_by_position(ds_cohort, size=window_size)
    ds_cohort_div["alt_hom"] = (("windows", "cohorts"), non_ref_pop(ds_cohort_div)[:,:,0])
    chrom_df = pd.DataFrame({"window_start": list(range(0, len(ds_cohort_div.window_start)*window_size, window_size))})
    chrom_df["divergence"] = ds_cohort_div.alt_hom/len(ds["samples"])
    ds_cohort_pi = ds_cohort.sel(variants = ((ds_cohort.cohort_allele_count[:,0] >= 1).sum(axis=1) >= 2).compute())
    if len(ds_cohort_pi.variants) < 10:
        chrom_df["pi"] = np.nan
    else:
        ds_cohort_pi = sg.window_by_position(ds_cohort_pi, size=window_size) 
        ds_cohort_pi = (sg.diversity(ds_cohort_pi))
        chrom_df["pi"] = pd.Series(ds_cohort_pi.stat_diversity[:,0])
    chrom_df["chrom"] = c
    df_l.append(chrom_df)
df_pi = pd.concat(df_l)
intervals_callable = pos_windows(bed_files, window_size, bed_files["chrom"].unique())
output_df = pd.merge(df_pi, intervals_callable, on=["chrom", "window_start"])
output_df["chr_type"] = output_df["chrom"].map(dict(zip(regions_df.CONTIG_ID, regions_df.chr_type)))
output_df["species"] = long_form
output_df.to_csv("{}{}_{}kb_pi.txt".format(args.o, long_form, args.w), sep="\t", index=False)
