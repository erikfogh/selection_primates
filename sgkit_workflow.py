from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

starting_dir = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/{}/filteredVCF/bcf_step1/"
starting_dir_gvcf = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/{}/gVCF/"
metadata_path = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
zarr_dir = "zarr_data/"

def generate_zarr(focus_chr, all_chr, out_path):
    inputs = [all_chr]
    o = out_path+focus_chr+"/call_genotype"
    outputs = [o]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "120:00:00",
        "account": "baboondiversity"
    }
    spec = """
    bcftools view -Ou -r {focus_chr} {all_chr} > {temp_bcf}
    bcftools index {temp_bcf}
    vcf2zarr convert {temp_bcf} {out_path}
    rm {temp_bcf}
    rm {temp_bcf}.csi
    """.format(focus_chr=focus_chr, all_chr=all_chr,
               temp_bcf=out_path+focus_chr+"_temp.bcf", out_path=out_path+focus_chr)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_vcf(idx, target):
    #print(target.spec)
    filename = target.spec.split("/")[-2].split(".")[0]+"_"+target.spec.split("/")[-1].split(".")[0]
    return 'vcf_zarr_{}'.format(filename)


# This section goes through all metadata folders and identifies the corresponding VCFs to use.

metadata_folders = glob.glob(metadata_path+"*_individuals.txt")
#print(metadata_folders)

for folder in metadata_folders:
    metadata_df = pd.read_csv(folder, sep="\t")
    short_form = folder.split("/")[-1].split("_")[0]
    regions_df = pd.read_csv(metadata_path+"{}_regions_and_batches.txt".format(short_form), sep="\t")
    for GVCF_FOLDER in metadata_df.GVCF_FOLDER.unique():
        #print(GVCF_FOLDER)
        reference = metadata_df.loc[metadata_df.GVCF_FOLDER == GVCF_FOLDER].REFERENCE_FOLDER.unique()
        all_chr_file = starting_dir.format(GVCF_FOLDER)+"{}_all_chr.sorted.bcf".format(GVCF_FOLDER)
        if not os.path.exists(all_chr_file):
            print("Filtered VCF does not exist ", GVCF_FOLDER)
            continue
        zarr_input = []
        # Chromosome X
        zarr_input.extend(regions_df.loc[(regions_df.END >= 1000000) &
                                      (regions_df.FEMALE_PLOIDY == 2) &
                                      (regions_df.MALE_PLOIDY == 1) &
              (regions_df.REFERENCE_FOLDER == reference[0])].CONTIG_ID.unique()[:2])
        # Autosomes
        zarr_input.extend(regions_df.loc[(regions_df.END >= 1000000) &
                                      (regions_df.FEMALE_PLOIDY == 2) &
                                      (regions_df.MALE_PLOIDY == 2) &
              (regions_df.REFERENCE_FOLDER == reference[0])].CONTIG_ID.unique()[:2])
        out_path = zarr_dir+GVCF_FOLDER+"/"
        os.makedirs(out_path, exist_ok=True)
        gwf.map(generate_zarr, zarr_input, name=get_ID_vcf,
                extra={"all_chr": all_chr_file, "out_path": out_path})
        # # Chromosome Y needs a different input file as it isn't part of the filteredVCF input.
        # I'm shuttering it for now - need at the very least bed files and PAR identification to be useful
        # y_inputs = regions_df.loc[(regions_df.END >= 1000000) &
        #                               (regions_df.FEMALE_PLOIDY == 0) &
        #                               (regions_df.MALE_PLOIDY == 1) &
        #       (regions_df.REFERENCE_FOLDER == reference[0])]
        # zarr_input = y_inputs.CONTIG_ID.unique()[:2]
        # if len(y_inputs["BATCH"]) != 1:
        #     # if len(y_inputs["BATCH"]) > 1:
        #     #     print("Y chromosome is present in", len(y_inputs["BATCH"]), "batches")
        #     continue
        # batch = y_inputs["BATCH"].iloc[0]
        # y_file = starting_dir_gvcf.format(GVCF_FOLDER)+"{}_batch_{}_fploidy_0_mploidy_1_gt.gvcf.gz".format(GVCF_FOLDER, batch)
        # if not os.path.exists(y_file):
        #     print("Batch does not exist ", y_file)
        #     continue
        # gwf.map(generate_zarr, zarr_input, name=get_ID_vcf,
        #         extra={"all_chr": y_file, "out_path": out_path})


# test_dir = [{"refname": "Papio_cynocephalus_ssp", "out_path": "data/Papio_cynocephalus_ssp/"},
#             {"refname": "Gorilla_gorilla_ssp", "out_path": "data/Gorilla_gorilla_ssp/"}]

# gwf.map(generate_zarr, test_dir, name=get_ID_vcf)

# # GVCF implementation. Will use it to bash together the different anubis/mulatta aligned species.

# def generate_zarr_gvcf(focus_chr, chr_list, samples, out_path):
#     inputs = [chr_list]
#     o = out_path+focus_chr+"/call_genotype"
#     outputs = [o]
#     options = {
#         "cores": 4,
#         "memory": "30g",
#         "walltime": "120:00:00",
#         "account": "baboondiversity"
#     }
#     spec = """
#     bcftools merge -Ou -r {focus_chr} {chr_list} |
#     bcftools view -Ou -s {samples} |
#     bcftools view -Ou -m2 -M2 |
#     bcftools view -Ou -i 'MAF[0]>=0.1' |
#     bcftools +prune -Ob -m LD=0.5 -w 1000bp --random-seed 10 -o {temp_bcf}
#     bcftools index {temp_bcf}
#     vcf2zarr convert --force {temp_bcf} {out_path}
#     rm {temp_bcf}
#     rm {temp_bcf}.csi
#     """.format(focus_chr=focus_chr, chr_list=" ".join(chr_list), samples=",".join(samples),
#                temp_bcf=out_path+focus_chr+"2_temp.bcf",
#                out_path=out_path+focus_chr)
#     # print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def get_ID_vcf_gvcf(idx, target):
#     #print(target.spec)
#     filename = target.spec.split("/")[-2].split(".")[0]+"_"+target.spec.split("/")[-1].split(".")[0]
#     return 'vcf_zarr_{}'.format(filename)


# metadata_folders_pca = ['/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/Papio_individuals.txt',
#                         '/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/Macaca_individuals.txt'] # glob.glob(metadata_path+"Pa*_individuals.txt")
# starting_dir_pca = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/{}/gVCF/"


# for folder in metadata_folders_pca:
#     metadata_df = pd.read_csv(folder, sep="\t")
#     short_form = folder.split("/")[-1].split("_")[0]
#     regions_df = pd.read_csv(metadata_path+"{}_regions_and_batches.txt".format(short_form), sep="\t")
#     metadata_df = metadata_df.loc[metadata_df.REFERENCE_FOLDER.isin(["Papio_anubis_ssp", "Macaca_mulatta_ssp"])]
#     for REFERENCE_FOLDER in metadata_df.REFERENCE_FOLDER.unique():
#         reference = REFERENCE_FOLDER
#         gvcfs = metadata_df.GVCF_FOLDER.unique()
#         print(reference)
#         zarr_input = []
#         regions_df = regions_df.loc[regions_df.REFERENCE_FOLDER == reference]
#         aut_and_x = regions_df.loc[(regions_df.FEMALE_PLOIDY == 2) &
#                                    (regions_df.END >= 1000000)]
#         inds = metadata_df.loc[(metadata_df.AVG_COVERAGE_A >= 20)].GVCF_ID
#         print(len(inds))
#         # print(aut_and_x)
#         for c in aut_and_x.CONTIG_ID.unique():
#             # Batches are either aut/Y-linked or X-linked.
#             mploidy = aut_and_x.loc[aut_and_x.CONTIG_ID == c].MALE_PLOIDY.iloc[0]
#             b = aut_and_x.loc[aut_and_x.CONTIG_ID == c].BATCH.iloc[0]
#             #print(c, mploidy, b)
#             in_paths = [starting_dir_pca.format(x)+"{}_batch_{}_fploidy_2_mploidy_{}_gt.gvcf.gz".format(x, b, mploidy) for x in gvcfs]
#             d = {}
#             d["focus_chr"] = c
#             d["chr_list"] = in_paths
#             zarr_input.append(d)
#         out_path = zarr_dir+reference+"_pca/"
#         os.makedirs(out_path, exist_ok=True)
#         gwf.map(generate_zarr_gvcf, zarr_input, name=get_ID_vcf_gvcf,
#                 extra={"samples": inds, "out_path": out_path})

