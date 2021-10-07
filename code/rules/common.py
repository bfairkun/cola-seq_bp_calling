import pandas as pd
import os
import hashlib

###### Config file and sample sheets #####

samples = pd.read_csv(".test/samples.tsv" ,sep='\t', index_col=0, comment='#')
test_samples = pd.read_csv(".test/samples.tsv" ,sep='\t', index_col=0, comment='#')
# samples = pd.read_csv(config["samples"],sep='\t', index_col=0)

# Check if samples file is the included test samples.tsv file. The fastq
# samples in that file come from a bam that yi shared with me, so they need to
# be processed sort of uniquely (eg, no adapter trimming, since adapters have
# already been trimmed)
Are_samples_test_samples = (hashlib.md5(open(config["samples"], 'rb').read()).hexdigest() == '3fd2bc9a2c44229a5488a0cdd234cc3d')
# if Are_samples_test_samples:
#     print("Using the test samples file (.test/samples.tsv) which has an appropriate checksum")
# else:
#     print("Using")

BP_groups = samples['BP_group'].unique()

def require_at_least_one(filelist):
    existing = [file for file in filelist if os.path.isfile(file)]
    return existing if len(existing) else "non_existing_file"

def flatten_list(_2d_list):
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

def GetFlattenedBamAndFastqList(_pd_samplestsv_datafame):
    fastq_and_bams = _pd_samplestsv_datafame[['Bam_file', 'Fastq_R1_list', 'Fastq_R2_list']].to_numpy().flatten().tolist()
    return flatten_list([i.split(',') for i in fastq_and_bams if pd.isnull(i) == False])

def IsSampleFromTestConfig(_Sample, test_samples = test_samples):
    """
    Samples from test config were generated from bam files that have already
    been deduplicated by UMIs, with UMIs and adapters sequences removed from
    the read. These samples need to be processed differently. This helper
    function returns True or False based on if the Sample's fastq or bam
    filepath in the samples.tsv file matches any sample in the test samples.tsv
    file.
    """
    Blacklist_bam_and_fastq = GetFlattenedBamAndFastqList(test_samples)
    Sample_bam_and_fastq = GetFlattenedBamAndFastqList(samples.loc[_Sample])
    if len(set(Blacklist_bam_and_fastq) & set(Sample_bam_and_fastq)) > 0:
        return True
    else:
        return False

def GetBamList_For_BP_group(wildcards):
    return expand("Alignments/Final/{sample}.bam", sample = samples.loc[samples['BP_group'] == wildcards.BP_group].index.unique())

def GetBaiList_For_BP_group(wildcards):
    return expand("Alignments/Final/{sample}.bam.bai", sample = samples.loc[samples['BP_group'] == wildcards.BP_group].index.unique())


def GetFastqList_For_Sample(Read):
    """
    Return an input function that returns fastq files, either R1 or R2, based
    on Read argument which should match a column name in the samples.tsv file
    """
    def F(wildcards):
        return samples.loc[wildcards.sample][Read].split(',')
    return F

def IsSampleFromTestConfig_smk_input_function(wildcards):
    """
    return 'true' or 'false', the bash boolean
    """
    return str(IsSampleFromTestConfig(wildcards.sample)).lower()

def Get_fastq_params(wildcards):
    """
    In input function to get fastp extra parameters, optionally utilizing the
    'sample' wildcard. In this example, I check if the fastq is from the test
    fastq, which needs different fastp params since reads are already trimmed
    of adapters and UMIs in the test data
    """
    if IsSampleFromTestConfig(wildcards.sample) == True:
        return "-L -A -G"
    # elif wildcards.sample == "Example"; return SomeOtherString
    else:
        return "--umi --umi_loc read2 --umi_len 6"

