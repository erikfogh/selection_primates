from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
zarr_path = "human_ref_zarr_20x_inds/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
feature_bed_path = "data/hs1_CATGenes.bed"
out_path = "results/human_window_stats_20x_inds/"
window_size = 25 # In kb
interval_size = 100 # Also in kb
fst_cutoff = 0.1

def analyse_zarr_pi(zarr_dir, window_size, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-2]
    outputs = ["{}{}_{}kb_pi.txt".format(out_path, long_form, window_size)]
    options = {
        "cores": 4,
        "memory": "28g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/human_zarr_analysis_pi.py -i {zarr_dir} -w {window_size} -o {out_path}
    """.format(zarr_dir=zarr_dir, window_size=window_size,out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def annotate_df(out_path, feature_bed, window_size):
    inputs = [out_path]
    outputs = [out_path[:-4]+"_annotated.txt"]
    options = {
        "cores": 1,
        "memory": "8g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/human_zarr_annotate.py -i {out_path} -g {feature_bed} -w {window_size}
    """.format(out_path=out_path, feature_bed=feature_bed, window_size=window_size)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# This is an old version, but I might use it later for pi estimates based on sub-pops.
def analyse_zarr_fst_pi(zarr_dir, metadata, window_size, fst_cutoff, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-2]
    outputs = ["{}{}_{}kb_pi_fst{}.txt".format(out_path, long_form, window_size, fst_cutoff)]
    options = {
        "cores": 13,
        "memory": "132g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_fst_pi.py -i {zarr_dir} -m {metadata} -w {window_size} -c {fst_cutoff} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size,
               fst_cutoff=fst_cutoff, out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_analyse_pi(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_pi'.format(species)

def get_ID_annotate_pi(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[3].split(" ")[0]
    return '{}_annotate'.format(species)

# Do all analysis in one python script.

map_input = []
for x in glob.glob(zarr_path+"*")[:]:
    d = {}
    d["zarr_dir"] = x+"/zarr"
    d["window_size"] = window_size
    d["out_path"] = out_path
    map_input.append(d)
pi_dfs = gwf.map(analyse_zarr_pi, map_input, name=get_ID_analyse_pi)

gwf.map(annotate_df, pi_dfs.outputs, extra={"feature_bed": feature_bed_path,
                                            "window_size": 100},
                                            name=get_ID_annotate_pi)
