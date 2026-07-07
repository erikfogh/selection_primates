from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

# Paths and variables
chain_file = "data/hg38ToGCA_009914755.4.over.chain.gz"
chain_file_hs = "data/hg38ToHs1.over.chain.gz"
map_input = [{"in_bed": "data/recomb_bed.bed", "chain_file": chain_file, "out_bed": "data/lifted_recomb.bed"},
             {"in_bed": "data/phast_10kb.bed", "chain_file": chain_file, "out_bed": "data/lifted_phast.bed"},
             {"in_bed": "data/recomb_bed.bed", "chain_file": chain_file_hs, "out_bed": "data/hs1_lifted_recomb.bed"},
             {"in_bed": "data/phast_10kb.bed", "chain_file": chain_file_hs, "out_bed": "data/hs1_lifted_phast.bed"}]

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

gwf.map(human_bed_lift, map_input)

