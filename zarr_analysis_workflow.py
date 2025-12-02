from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
zarr_path = "zarr_20x_inds/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
out_path = "results/window_stats_20x_inds/"
window_size = 10 # In kb


def analyse_zarr_het_hom(zarr_dir, metadata, window_size, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-1]
    outputs = ["{}{}_{}kb_het_hom.txt".format(out_path, long_form, window_size)]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "60:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_het_hom.py -i {zarr_dir} -m {metadata} -w {window_size} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size, out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_analyse_het(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_het_hom'.format(species)

# Do all analysis in one python script.

map_input = []
for x in glob.glob(zarr_path+"*")[:]:
    d = {}
    d["zarr_dir"] = x
    d["metadata"] = metadata_path
    d["window_size"] = window_size
    d["out_path"] = out_path
    map_input.append(d)
gwf.map(analyse_zarr_het_hom, map_input, name=get_ID_analyse_het)
