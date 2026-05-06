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
table_desc = "~/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/plots/SupTable_Sample_Stats_wGT_QC.tsv"

size_cutoff = 1000000 # 1Mb
chr_cut_hetsep = 5
chr_cut_grid = 5
ind_cut = 1
decode_cut = 5


hetsep_dir = "steps/multihetsep"
cobraa_dir = "steps/cobraa"

# I do not distinguish between ind mask/call in this workflow.


def generate_multihetsep(refname, vcf, contig, mask, callability):
    out_path = contig+".txt"
    inputs = [vcf, mask, callability]
    outputs = ["steps/multihetsep/"+refname+"/"+out_path]
    options = {
        "cores": 1,
        "memory": "8g",
        "walltime": "4:00:00",
        "account": "primatediversity"
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


def cobraa_run(multihetsep_l, ts, te, chrom, mem=24):
    refname = multihetsep_l[0].split("/")[-2]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_D50_ts{}_te{}_".format(ts, te)
    inputs = multihetsep_l
    outputs = out_path+"final_parameters.txt"
    options = {
        "cores": 3,
        "memory": "{}g".format(mem),
        "walltime": "4:00:00",
        "account": "primatediversity"
    }
    spec = """
    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -ts {ts} -te {te} -its 20 -spread_1 0.025 -spread_2 50 
    """.format(infiles=" ".join(multihetsep_l), out_path=out_path, ts=ts, te=te)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def cobraa_run_unstructured(multihetsep_l, chrom, mem=24):
    refname = multihetsep_l[0].split("/")[-2]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_"
    inputs = multihetsep_l
    outputs = out_path+"final_parameters.txt"
    options = {
        "cores": 3,
        "memory": "{}g".format(mem),
        "walltime": "4:00:00",
        "account": "primatediversity"
    }
    spec = """
    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -its 20 -spread_1 0.025 -spread_2 50
    """.format(infiles=" ".join(multihetsep_l), out_path=out_path)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def cobraa_decode(param_l, multihetsep, chrom, mem=24):
    refname = multihetsep.split("/")[-2]
    chr_name = multihetsep.split("/")[-1][:-4]
    out_path = "steps/cobraa/"+refname+"/"+chrom+"_"+chr_name+"_decode.txt"
    inputs = param_l
    outputs = out_path
    joined_params = " ".join(param_l)
    options = {
        "cores": 2,
        "memory": "{}g".format(mem),
        "walltime": "4:00:00",
        "account": "primatediversity"
    }
    spec = """
    best_file_params=$(python scripts/cobraa_decode_pick.py -i {joined_params})
    python cobraa/cobraa.py -in {infiles} -o {out_path} -D 50 -b 25 -its 1 -thresh 1 -spread_1 0.025 -spread_2 50 \
    -decode -decode_downsample 400 -path -lambda_B_segments 50*1 \
    $best_file_params
    """.format(joined_params=joined_params,
               infiles=multihetsep, out_path=out_path)
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
    filename = target.spec.split("/")[-2]+target.spec.split("/")[-1].split("_")[0]+"_"+target.spec.split("/")[-1].split("_")[3]+target.spec.split("/")[-1].split("_")[4]
    #print(filename)
    return 'structured_cobraa_{}'.format(filename)

def get_ID_decode_cobraa(idx, target):
    #print(target.spec)
    filename = target.spec.split("-o")[-1].split("/")[-2]+"_"+target.outputs.split("/")[-1]
    #print(filename)
    return 'decode_cobraa_{}'.format(filename)


# Generating supporting directories
os.makedirs(hetsep_dir, exist_ok=True)
os.makedirs(cobraa_dir, exist_ok=True)


metadata_table = pd.read_csv(table_desc, sep="\t")
metadata_20x_filt = metadata_table.loc[(metadata_table.finalQC != "fail")
                              & (metadata_table.cov_chrA >= 20)
                              & (metadata_table.remove_as_relative != True)
                              & (metadata_table.remove_manual != True)
                              & (~metadata_table.ID.str.startswith("SAMEA11633"))
                             ]

skipped_cases = []
for g in metadata_20x_filt.genus.unique()[:3]:
    # Identify IDs
    species_metadata = metadata_20x_filt.loc[metadata_20x_filt.genus == g]
    species_metadata = species_metadata.loc[species_metadata.cov_chrX >= 5]
    # I only want females as they have a diploid chrX.
    female_df = species_metadata.loc[(species_metadata.gSEX == "F")].sort_values(by="cov_chrA", ascending=False)
    # # Restricted analysis to females - below is the way to include males if you want to run non-X analysis.
    # male_df = species_metadata.loc[~(species_metadata.ID.isin(female_df.ID))].sort_values(by="cov_chrA", ascending=False)
    # sorted_df = pd.concat([female_df, male_df])
    sorted_df = female_df
    if len(sorted_df) > 0:
        sorted_input = sorted_df
    else:
        print("No usable ind for", g)
        continue
    region_metadata = pd.read_csv(metadata_dir+g+"_regions_and_batches.txt",
                           sep="\t").sort_values(by="END")
    # Go through every unique genotype calling set.
    for gvcf_folder in sorted_input.species_genotyping.unique()[:]:
        # Pick the (currently one) best individuals in the metadata.
        picked_inds = sorted_input.loc[sorted_input.species_genotyping == gvcf_folder].ID.iloc[:ind_cut]
        # Choose all paths based on the regions file.
        batch_name = "{}/filteredVCF/all_samples/bcf_step1/{}_batch{}_fploidy2_mploidy{}.bcf"
        gvcfs_names = starting_dir+batch_name
        mask_names = starting_dir+"{}/filteredVCF/all_samples/pos_bed_cov_based/{}_batch{}_fploidy2_mploidy{}.bed"
        # Workflow selecting all autosomal chromosomes above 1Mb in size.
        aut_l = region_metadata.loc[(region_metadata.FEMALE_PLOIDY == 2) &
                                    (region_metadata.MALE_PLOIDY == 2) &
                                     (region_metadata.END >= size_cutoff)][:chr_cut_hetsep]
        x_l = region_metadata.loc[(region_metadata.FEMALE_PLOIDY == 2) &
                                  (region_metadata.MALE_PLOIDY == 1) &
                                     (region_metadata.END >= size_cutoff)][:chr_cut_hetsep]
        aut_and_x = pd.concat([aut_l, x_l])
        # print(aut_and_x)
        for ind in picked_inds:
            multihetsep_args = []
            for c in aut_and_x.CONTIG_ID.unique():
            # Batches are either aut/Y-linked or X-linked.
                mploidy = aut_and_x.loc[aut_and_x.CONTIG_ID == c].MALE_PLOIDY.iloc[0]
                b = aut_and_x.loc[aut_and_x.CONTIG_ID == c].BATCH.iloc[0]
                if os.path.exists(gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy)) and os.path.exists(mask_names.format(gvcf_folder, gvcf_folder, b, mploidy)):
                    os.makedirs(hetsep_dir+"/"+ind, exist_ok=True)
                    ind_dir = {}
                    ind_dir["refname"] = ind
                    ind_dir["vcf"] = gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    ind_dir["contig"] = c
                    ind_dir["mask"] = mask_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    ind_dir["callability"] = mask_names.format(gvcf_folder, gvcf_folder, b, mploidy)
                    multihetsep_args.append(ind_dir)
                else:
                    skipped_cases.append(gvcf_folder)
                    # print(gvcfs_names.format(gvcf_folder, gvcf_folder, b, mploidy), "does not exist")
                    continue
            multihetsep_o = gwf.map(generate_multihetsep, multihetsep_args, name=get_ID_vcf)
            hetsep_l, hetseps_x_l = [], []
            os.makedirs(cobraa_dir+"/"+ind, exist_ok=True)
            # Aut run
            c_list = []
            for c in aut_l.CONTIG_ID.unique()[:chr_cut_grid]:
                c_list.append("steps/multihetsep/{}/{}.txt".format(ind, c))
                # hetsep_l.append({"multihetsep_l": ["steps/multihetsep/{}/{}.txt".format(i, c)],
                # "chrom": c})
            hetsep_l.append({"multihetsep_l": c_list, "chrom": "aut_PSMC",
                             "mem": 24})
            x_list = []
            for c in x_l.CONTIG_ID.unique()[:chr_cut_grid]:
                x_list.append("steps/multihetsep/{}/{}.txt".format(ind, c))
            hetseps_x_l.append({"multihetsep_l": x_list, "chrom": "chrX_PSMC",
                             "mem": 16})
            # Unstructured aut
            gwf.map(cobraa_run_unstructured, hetsep_l,
                    name=get_ID_unstructured_cobraa)
            # Structured aut search
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
                # Decode can only take one file, so iterate through the various contigs
                decode_l = []
                for c in aut_l.CONTIG_ID.unique()[:decode_cut]:
                    decode_d = {}
                    decode_d["param_l"] = i_structured.outputs
                    decode_d["chrom"] = d["chrom"]
                    decode_d["multihetsep"] = "steps/multihetsep/{}/{}.txt".format(ind, c)
                    decode_d["mem"] = 24
                    decode_l.append(decode_d)
                
                gwf.map(cobraa_decode, decode_l,
                name=get_ID_decode_cobraa)

            # Unstructured X
            gwf.map(cobraa_run_unstructured, hetseps_x_l,
                    name=get_ID_unstructured_cobraa)
            # Structured X search
            for d in hetseps_x_l:
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
                # Decode can only take one file, so iterate through the various contigs
                decode_l = []
                for c in x_l.CONTIG_ID.unique()[:decode_cut]:
                    decode_d = {}
                    decode_d["param_l"] = i_structured.outputs
                    decode_d["chrom"] = d["chrom"]
                    decode_d["multihetsep"] = "steps/multihetsep/{}/{}.txt".format(ind, c)
                    decode_d["mem"] = 24
                    decode_l.append(decode_d)
                gwf.map(cobraa_decode, decode_l,
                name=get_ID_decode_cobraa)

print("Problems with", set(skipped_cases))

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
