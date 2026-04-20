from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
zarr_path = "zarr_20x_inds/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
chain_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/references49/alignment/chains/"
feature_bed_path = "data/hs1_CATGenes.bed"
out_path = "results/window_stats_20x_inds/"
out_path_lift = "results/lifted_window_stats_20x_inds/"
window_size = 10 # In kb

window_size_TD = 100
step_size_TD = 10

table_desc = "~/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/plots/SupTable_Sample_Stats_wGT_QC.tsv"
metadata_table = pd.read_csv(table_desc, sep="\t")

metadata_20x_filt = metadata_table.loc[(metadata_table.finalQC != "fail")
                              & (metadata_table.cov_chrA >= 20)
                              & (metadata_table.remove_as_relative != True)
                              & (metadata_table.remove_manual != True)
                              & (~metadata_table.ID.str.startswith("SAMEA11633"))
                             ]

count_sub = metadata_20x_filt.loc[~metadata_20x_filt.cov_chrX.isna()][["gSEX", "group", "species_genotyping", "species"]].value_counts().reset_index()
used_species = count_sub.loc[(count_sub.gSEX == "F") & (count_sub["count"] >= 3)].species_genotyping.unique()[:]
zarr_species = glob.glob(zarr_path+"*")

def analyse_zarr_pi(zarr_dir, metadata, window_size, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-2]
    outputs = ["{}{}_{}kb_pi.txt".format(out_path, long_form, window_size),
                "{}{}_{}kb_pi.bed".format(out_path, long_form, window_size)]
    options = {
        "cores": 4,
        "memory": "28g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_pi.py -i {zarr_dir} -m {metadata} -w {window_size} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size,out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def analyse_zarr_divergence_pi(zarr_dir, metadata, window_size, fst_cutoff, out_path):
#     inputs = [zarr_dir]
#     long_form = zarr_dir.split("/")[-2]
#     outputs = ["{}{}_biggest_pop.txt".format(out_path, long_form),
#                "{}{}_{}kb_pi.txt".format(out_path, long_form, window_size)]
#     options = {
#         "cores": 4,
#         "memory": "28g",
#         "walltime": "12:00:00",
#         "account": "baboondiversity"
#     }
#     spec = """
#     python scripts/zarr_analysis_pi.py -i {zarr_dir} -m {metadata} -w {window_size} -o {out_path}
#     """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size,out_path=out_path)
#     # print(outputs, spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def analyse_zarr_Tajima_D(zarr_dir, metadata, window_size, step_size, out_path):
    inputs = [zarr_dir]
    long_form = zarr_dir.split("/")[-2]
    outputs = ["{}{}_{}kb_{}step_Tajima_D.bed".format(out_path, long_form, window_size, step_size)]
    options = {
        "cores": 4,
        "memory": "28g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/zarr_analysis_tajimaD.py -i {zarr_dir} -m {metadata} -w {window_size} -s {step_size} -o {out_path}
    """.format(zarr_dir=zarr_dir, metadata=metadata, window_size=window_size, step_size=step_size,out_path=out_path)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def human_bed_lift(in_bed, chain_file, out_bed):
    inputs = in_bed
    outputs = out_bed
    options = {
        "cores": 4,
        "memory": "28g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    CrossMap bed {chain_file} {in_bed} {out_bed}
    """.format(in_bed=in_bed, chain_file=chain_file, out_bed=out_bed)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bed_annotate_window_df(out_path, feature_bed, window_size):
    inputs = [out_path]
    outputs = out_path[:-4]+"_windowed_annotated.txt"
    options = {
        "cores": 2,
        "memory": "16g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/bed_annotate_window.py -i {out_path} -g {feature_bed} -w {window_size}
    """.format(out_path=out_path, feature_bed=feature_bed, window_size=window_size)
    # print(outputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_analyse_pi(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_pi'.format(species)


def get_ID_analyse_TD(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[2].split(" ")[0]
    return '{}_Tajima'.format(species)

def get_ID_analyse_bed(idx, target):
    #print(target.spec)
    out_bed = target.outputs.split("/")[-1]
    return '{}'.format(out_bed)

def get_ID_annotate_pi(idx, target):
    #print(target.spec)
    species = target.spec.split("/")[3].split(" ")[0]
    return '{}_annotate'.format(species)

# Mapping over the different gwf workflows

# map_input = []
# for x in zarr_species:
#     d = {}
#     d["zarr_dir"] = x+"/zarr"
#     d["metadata"] = metadata_path
#     d["window_size"] = window_size
#     d["out_path"] = out_path
#     map_input.append(d)
# gwf.map(analyse_zarr_pi, map_input, name=get_ID_analyse_pi)

# map_input = []
# for x in zarr_species:
#     s = x.split("/")[-1]
#     chain_ref = metadata_20x_filt.loc[metadata_20x_filt.species_genotyping == s].reference.iloc[0]
#     d = {}
#     d["in_bed"] = "{}{}_{}kb_pi.bed".format(out_path, s, window_size)
#     d["chain_file"] = chain_path+chain_ref+"_vs_Homo_sapiens.chain.gz"
#     d["out_bed"] = "{}{}_{}kb_pi.bed".format(out_path_lift, s, window_size)
#     map_input.append(d)
# bed_lift = gwf.map(human_bed_lift, map_input, name=get_ID_analyse_bed)

# gwf.map(bed_annotate_window_df, bed_lift.outputs, extra={"feature_bed": feature_bed_path,
#                                             "window_size": 100},
#                                             name=get_ID_annotate_pi)



map_input = []
for s in used_species:
    d = {}
    d["zarr_dir"] = zarr_path+s+"/zarr"
    d["metadata"] = metadata_path
    d["window_size"] = window_size_TD
    d["step_size"] = step_size_TD
    d["out_path"] = out_path
    map_input.append(d)
gwf.map(analyse_zarr_Tajima_D, map_input, name=get_ID_analyse_TD)

map_input = []
for s in used_species:
    chain_ref = metadata_20x_filt.loc[metadata_20x_filt.species_genotyping == s].reference.iloc[0]
    d = {}
    d["in_bed"] = "{}{}_{}kb_{}step_Tajima_D.bed".format(out_path, s, window_size_TD, step_size_TD)
    d["chain_file"] = chain_path+chain_ref+"_vs_Homo_sapiens.chain.gz"
    d["out_bed"] = "{}{}_{}kb_{}step_Tajima_D.bed".format(out_path_lift, s, window_size_TD, step_size_TD)
    map_input.append(d)
gwf.map(human_bed_lift, map_input, name=get_ID_analyse_bed)
