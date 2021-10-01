import pandas as pd
import os

###### Config file and sample sheets #####

samples = pd.read_csv(".test/samples.tsv" ,sep='\t', index_col=0, comment='#')
# samples = pd.read_csv(config["samples"],sep='\t', index_col=0)

BP_groups = samples['BP_group'].unique()


def GetBamList_For_BP_group(wildcards):
    return samples.loc[samples['BP_group'] == wildcards.BP_group]['Bam_file']

def GetBaiList_For_BP_group(wildcards):
    return [i + ".bai" for i in samples.loc[samples['BP_group'] == wildcards.BP_group]['Bam_file']]

