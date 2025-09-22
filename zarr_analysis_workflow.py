from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
zarr_path = "zarr_data/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
out_path = "results/window_stats/"
window_size = 100 # In kb


def analyse_zarr_het(zarr_dir, metadata, window_size, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-1]
    outputs = ["{}{}_{}kb_het.txt".format(out_path, long_form, window_size)]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "60:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_het.py -i {zarr_dir} -m {metadata} -w {window_size} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size, out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_analyse_het(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_het'.format(species)


# def analyse_zarr_quadratic(focus_chr, all_chr, out_path):
#     inputs = [all_chr]
#     o = out_path+focus_chr+"/call_genotype"
#     outputs = [o]
#     options = {
#         "cores": 4,
#         "memory": "30g",
#         "walltime": "120:00:00",
#         "account": "baboondiversity"
#     }
#     spec = """
#     python scripts/zarr_analysis.py -i {zarr_dir} -m {metadata} -n {ref_name} -w {window_size} -o {out_path}
#     """.format(focus_chr=focus_chr, all_chr=all_chr,
#                temp_bcf=out_path+focus_chr+"_temp.bcf", out_path=out_path+focus_chr)
#     # print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Do all analysis in one python script.

map_input = []
for x in glob.glob(zarr_path+"*")[10:20]:
    d = {}
    d["zarr_dir"] = x
    d["metadata"] = metadata_path
    d["window_size"] = window_size
    d["out_path"] = out_path
    map_input.append(d)
gwf.map(analyse_zarr_het, map_input, name=get_ID_analyse_het)