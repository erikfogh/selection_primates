from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
zarr_path = "zarr_20x_inds/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
out_path = "results/window_stats_20x_inds/"
window_size = 100 # In kb
fst_cutoff = 0.1


def analyse_zarr_fst_pi(zarr_dir, metadata, window_size, fst_cutoff, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-2]
    outputs = ["{}{}_{}kb_pi_fst{}.txt".format(out_path, long_form, window_size, fst_cutoff)]
    options = {
        "cores": 2,
        "memory": "16g",
        "walltime": "60:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_fst_pi.py -i {zarr_dir} -m {metadata} -w {window_size} -c {fst_cutoff} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size,
               fst_cutoff=fst_cutoff, out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_analyse_het(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_fst_pi'.format(species)

# Do all analysis in one python script.

map_input = []
for x in glob.glob(zarr_path+"*")[:5]:
    d = {}
    d["zarr_dir"] = x+"/zarr"
    d["metadata"] = metadata_path
    d["window_size"] = window_size
    d["fst_cutoff"] = fst_cutoff
    d["out_path"] = out_path
    map_input.append(d)
gwf.map(analyse_zarr_fst_pi, map_input, name=get_ID_analyse_het)
