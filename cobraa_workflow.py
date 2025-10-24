from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import numpy as np
import glob

gwf = Workflow()

# Paths. Autosome implementation. This setup is optimized for the primate_diversity directory setup.
# So it assumes a variety of things regarding file structure to make it sleeker.
# Multihetsep requires called sites, so I also filter for heterozygous sites as that is all it needs.
starting_dir = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024/"
metadata_dir = "/home/eriks/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/"
size_cutoff = 10000000


hetsep_dir = "steps/multihetsep"
cobraa_dir = "steps/cobraa"

# I do not distinguish between ind mask/call in this workflow.


def generate_multihetsep(refname, vcf, contig, mask, callability):
    out_path = contig+".txt"
    inputs = [vcf, mask, callability]
    outputs = ["steps/multihetsep/"+refname+"/"+out_path]
    options = {
        "cores": 2,
        "memory": "14g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    bcftools view -Ou -s {refname} -r {contig} {vcf} | bcftools view -Oz -v snps -g ^miss > steps/multihetsep/{refname}/{contig}_temp.vcf.gz
    grep "{contig}" {mask} > steps/multihetsep/{refname}/{contig}_temp.bed
    python cobraa/msmc-tools/generate_multihetsep.py --mask=steps/multihetsep/{refname}/{contig}_temp.bed steps/multihetsep/{refname}/{contig}_temp.vcf.gz>steps/multihetsep/{refname}/{out_path}
    rm steps/multihetsep/{refname}/{contig}_temp.vcf.gz
    rm steps/multihetsep/{refname}/{contig}_temp.bed
    """.format(vcf=vcf, refname=refname, contig=contig, mask=mask, out_path=out_path)
    #print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def cobraa_run(multihetsep_l, ts, te, chrom, mem=30):
    refname = multihetsep_l[0].split("/")[-2]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_D50_ts{}_te{}_".format(ts, te)
    inputs = multihetsep_l
    outputs = out_path+"final_parameters.txt"
    options = {
        "cores": 3,
        "memory": "{}g".format(mem),
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -ts {ts} -te {te} -its 20 -spread_1 0.025 -spread_2 50 
    """.format(infiles=" ".join(multihetsep_l), out_path=out_path, ts=ts, te=te)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def cobraa_run_unstructured(multihetsep_l, chrom, mem=30):
    refname = multihetsep_l[0].split("/")[-2]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_"
    inputs = multihetsep_l
    outputs = out_path+"final_parameters.txt"
    options = {
        "cores": 3,
        "memory": "{}g".format(mem),
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -its 20 -spread_1 0.025 -spread_2 50 
    """.format(infiles=" ".join(multihetsep_l), out_path=out_path)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def cobraa_decode(multihetsep_l, chrom, mem=120):
    refname = multihetsep_l[0].split("/")[-2]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_"
    inputs = multihetsep_l
    outputs = [out_path+"final_parameters.txt"]
    options = {
        "cores": 6,
        "memory": "{}g".format(mem),
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """

    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -its 20 -spread_1 0.025 -spread_2 50 
    """.format(infiles=" ".join(multihetsep_l), out_path=out_path)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_vcf(idx, target):
    #print(target.spec)
    filename = target.spec.split("/")[-2]+"_"+target.spec.split("/")[-1].split("_temp")[0]
    return 'vcf_transform_{}'.format(filename)

def get_ID_unstructured_cobraa(idx, target):
    #print(target.spec)
    filename = target.spec.split("/")[-2]+"_"+target.spec.split("/")[-1].split("_")[0]
    return 'unstructured_cobraa_{}'.format(filename)

def get_ID_structured_cobraa(idx, target):
    #print(target.spec)
    filename = target.spec.split("/")[-1].split("_")[0]+"_"+target.spec.split("/")[-1].split("_")[2]+target.spec.split("/")[-1].split("_")[3]
    #print(filename)
    return 'structured_cobraa_{}'.format(filename)


# Generating supporting directories
os.makedirs(hetsep_dir, exist_ok=True)
os.makedirs(cobraa_dir, exist_ok=True)

metadata_dirs = glob.glob(metadata_dir+"*_individuals.txt")

chr_cut = 2

for d in metadata_dirs[:1]:
    # Identify IDs
    dir_metadata = pd.read_csv(d, sep="\t")
    dir_metadata["gss"] = dir_metadata.GENUS+"_"+dir_metadata.SPECIES+"_"+dir_metadata.SUBSPECIES
    # Slightly hacky way of preferring females and those with high coverage.
    short_species = d.split("/")[-1].split("_")[0]
    female_df = dir_metadata[pd.to_numeric(dir_metadata['AVG_COVERAGE_X'], errors='coerce').notnull()]
    female_df = female_df.loc[(female_df.GENETIC_SEX == "F") & (female_df.AVG_COVERAGE_A >= 10)].sort_values(by="AVG_COVERAGE_A", ascending=False)
    male_df = dir_metadata[pd.to_numeric(dir_metadata['AVG_COVERAGE_A'], errors='coerce').notnull()]
    male_df = male_df.loc[~(male_df.GVCF_ID.isin(female_df.GVCF_ID)) & (male_df.AVG_COVERAGE_A >= 10)].sort_values(by="AVG_COVERAGE_A", ascending=False)
    sorted_df = pd.concat([female_df, male_df])
    if len(sorted_df) > 0:
        sorted_input = sorted_df
    else:
        print("No usable ind for", short_species)
        continue
    region_metadata = pd.read_csv(metadata_dir+short_species+"_regions_and_batches.txt",
                           sep="\t")
    # Go through every unique genotype calling set.
    for gvcf_folder in sorted_input.GVCF_FOLDER.unique()[:1]:
        print("Using", gvcf_folder)
        # Pick the (currently one) best individuals in the metadata.
        picked_inds = sorted_input.loc[sorted_input.GVCF_FOLDER == gvcf_folder].GVCF_ID.iloc[:1]
        # Choose all paths based on the regions file.
        batch_name = "{}/filteredVCF/bcf_step1/{}_batch{}_fploidy2_mploidy{}.bcf"
        gvcfs_names = starting_dir+batch_name
        mask_names = starting_dir+"{}/filteredVCF/pos_bed_cov_based/{}_batch{}_fploidy2_mploidy{}.bed"
        multihetsep_args = []
        # Workflow selecting all autosomal chromosomes above 1Mb in size.
        aut_l = region_metadata.loc[(region_metadata.FEMALE_PLOIDY == 2) &
                                    (region_metadata.MALE_PLOIDY == 2) &
                                     (region_metadata.END >= size_cutoff)][:chr_cut]
        x_l = region_metadata.loc[(region_metadata.FEMALE_PLOIDY == 2) &
                                  (region_metadata.MALE_PLOIDY == 1) &
                                     (region_metadata.END >= size_cutoff)][:chr_cut]
        aut_and_x = pd.concat([aut_l, x_l])
        # print(aut_and_x)
        for c in aut_and_x.CONTIG_ID.unique():
            # Batches are either aut/Y-linked or X-linked.
            mploidy = aut_and_x.loc[aut_and_x.CONTIG_ID == c].MALE_PLOIDY.iloc[0]
            b = aut_and_x.loc[aut_and_x.CONTIG_ID == c].BATCH.iloc[0]
            for i in picked_inds:
                if os.path.exists(gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy)):
                    os.makedirs(hetsep_dir+"/"+i, exist_ok=True)
                    ind_dir = {}
                    ind_dir["refname"] = i
                    ind_dir["vcf"] = gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    ind_dir["contig"] = c
                    ind_dir["mask"] = mask_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    ind_dir["callability"] = mask_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    multihetsep_args.append(ind_dir)
                else:
                    print(gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy), "does not exist")
        multihetsep_o = gwf.map(generate_multihetsep, multihetsep_args, name=get_ID_vcf)
        hetsep_l = []
        for i in picked_inds:
            os.makedirs(cobraa_dir+"/"+i, exist_ok=True)
            # Aut run
            c_list = []
            for c in aut_l.CONTIG_ID.unique():
                c_list.append("steps/multihetsep/{}/{}.txt".format(i, c))
                # hetsep_l.append({"multihetsep_l": ["steps/multihetsep/{}/{}.txt".format(i, c)],
                # "chrom": c})
            hetsep_l.append({"multihetsep_l": c_list, "chrom": "aut",
                             "mem": 60})

        #     x_list = []
        #     for b in region_metadata.loc[(region_metadata.FEMALE_PLOIDY == 2) &
        #                                  (region_metadata.MALE_PLOIDY == 1) &
        #                                  (region_metadata.END >= 1000000)].BATCH.unique():
        #         x_list.append("steps/multihetsep/{}/{}_batch{}.txt".format(i, i, b))
        #     hetsep_l.append({"multihetsep_l": x_list, "bnames": "_chrX"})
        # Unstructured
        gwf.map(cobraa_run_unstructured, hetsep_l,
                name=get_ID_unstructured_cobraa)
        for d in hetsep_l:
            i_l = []
            for i in range(10, 42, 6):
                for j in range(4, i-4, 6):
                    l = d.copy()
                    l["ts"] = j
                    l["te"] = i
                    i_l.append(l)
            i_structured = gwf.map(cobraa_run, i_l,
                                   name=get_ID_structured_cobraa)
            # Pick best and perform cobraa-path.
            print(i_structured.outputs)
            # gwf.map(cobraa_decode, i_structured.outputs)




# # Test of the old baboondiversity to compare.


# def generate_multihetsep_baboon(refname, vcf, mask, callability):
#     refname_batch = refname+"_"+vcf.split("/")[-1].split(".")[0]
#     multihetsep_name = refname_batch+".txt"
#     out_path = multihetsep_name
#     inputs = [vcf, mask, callability]
#     outputs = ["steps/multihetsep/"+refname+"/"+out_path]
#     options = {
#         "cores": 2,
#         "memory": "14g",
#         "walltime": "12:00:00",
#         "account": "baboondiversity"
#     }
#     spec = """
#     bcftools view -Ou -s {refname} {vcf} | bcftools view -Oz -g ^miss > steps/multihetsep/{refname}/{refname_batch}_temp.vcf.gz
#     python cobraa/msmc-tools/generate_multihetsep.py --negative_mask={mask} steps/multihetsep/{refname}/{refname_batch}_temp.vcf.gz>steps/multihetsep/{refname}/{out_path}
#     rm steps/multihetsep/{refname}/{refname_batch}_temp.vcf.gz
#     """.format(vcf=vcf, refname=refname, refname_batch=refname_batch, mask=mask, out_path=out_path)
#     #print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def cobraa_run_unstructured_baboon(multihetsep_l, bnames, mem=30):
#     refname = multihetsep_l[0].split("/")[-2]
#     out_path = "steps/cobraa/"+refname+bnames+"_"
#     inputs = multihetsep_l
#     outputs = [out_path+"final_parameters.txt"]
#     options = {
#         "cores": 3,
#         "memory": "{}g".format(mem),
#         "walltime": "12:00:00",
#         "account": "baboondiversity"
#     }
#     spec = """
#     python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 100 -its 20
#     """.format(infiles=" ".join(multihetsep_l), out_path=out_path)
#     #print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def get_ID_unstructured_cobraa_baboon(idx, target):
#     print(target.spec)
#     filename = target.spec.split("/")[-1].split("_")[0]+"_"+target.spec.split("/")[-1].split(" ")[0]
#     return 'unstructured_cobraa_{}'.format(filename)


# # for d in ["/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/data/Papio_metadata_with_clustering_sci.txt"]:
# #     # Identify IDs
# #     dir_metadata = pd.read_csv(d, sep =" ")
# #     picked_inds = ["PD_0214", "PD_0693"]
# #     for gvcf_folder in ["pass"]:
# #         vcf_name = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.phased.rehead.vcf.gz"
# #         mask_names = "/home/eriks/baboondiversity/data/callability_panu3_26_04_2021/chr{}sorted.bed.gz"
# #         multihetsep_args = []
# #         for b in ["{}".format(x) if x<= 20 else "X" for x in range(1, 22)]:
# #             for i in picked_inds:
# #                 os.makedirs(hetsep_dir+"/"+i, exist_ok=True)
# #                 ind_dir = {}
# #                 ind_dir["refname"] = i
# #                 ind_dir["vcf"] = vcf_name.format(b, b)
# #                 ind_dir["mask"] = mask_names.format(b)
# #                 ind_dir["callability"] = mask_names.format(b)
# #                 multihetsep_args.append(ind_dir)
# #         multihetsep_o = gwf.map(generate_multihetsep_baboon, multihetsep_args)
# #         hetsep_l = []
# #         for i in picked_inds:
# #             # Aut run
# #             b_list = []
# #             for b in range(1, 21):
# #                 b_list.append("steps/multihetsep/{}/{}_chr{}.txt".format(i, i, b))
# #                 hetsep_l.append({"multihetsep_l": ["steps/multihetsep/{}/{}_chr{}.txt".format(i, i, b)],
# #                 "bnames": "_chr{}".format(b)})
# #             hetsep_l.append({"multihetsep_l": b_list, "bnames": "_aut",
# #                              "mem": 60})

# #             x_list = []
# #             for b in ["X"]:
# #                 x_list.append("steps/multihetsep/{}/{}_chr{}.txt".format(i, i, b))
# #             hetsep_l.append({"multihetsep_l": ["steps/multihetsep/{}/{}_chr{}.txt".format(i, i, b)],
# #                 "bnames": "_chr{}".format(b)})
# #         # Unstructured
# #         gwf.map(cobraa_run_unstructured_baboon, hetsep_l,
# #                 name=get_ID_unstructured_cobraa_baboon)
