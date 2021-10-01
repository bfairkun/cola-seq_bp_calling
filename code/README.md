# Snakemake workflow: CoLa-seq_BP_calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/cola-seq_bp_calling.svg?branch=master)](https://travis-ci.org/snakemake-workflows/cola-seq_bp_calling)


## Authors

* Benjamin Fair (@bfairkun)

## Overview

The code is implemented as a snakemake. Some parameters for BP calling are set in `config/config/yaml`. Samples can be defined in `config/samples.tsv`. I will try to include rules for read mapping (from fastq files), or if you have already mapped reads, you can just define bam files in `config/samples.tsv`. Some rules create files that I have already included in the repo so that you probably don't need to redo those computations (eg, calculating a position weight matrices for published BPs) unless you want to modify those scripts/rules yourself for some reason, or, if for example, if you want to work with a different reference genome, then you may also need to rerun those rules. The version of `config/samples.tsv` as included in this repository will run the pipeline on test data that is also included in this repository, so the snakemake is already configured to run as is on the test data. The test data consists of all the reads from all libraries in Yi's manuscript that are around the URB1 gene (part of which is visualized in Fig1D in the manuscript)

In addition to BP mapping, I have included rules to tabulate BP usage at called BPs in each sample. I recommend grouping all samples together for BP calling to maximize BP discovery power, then tabulate sample specific BP counts. If you are working with vastly different tissue types, then you may consider grouping samples by tissue for BP calling.

## Usage

### Step 1: Install workflow and dependencies

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/cola-seq_bp_calling).

    git clone git@github.com:bfairkun/cola-seq_bp_calling.git

If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

Install snakemake and the workflow's other dependencies via conda/mamba. If conda/mamba isn't already installed, I recommend [installing miniconda](https://docs.conda.io/en/latest/miniconda.html) and then [install mamba](https://github.com/mamba-org/mamba) in your base environment. Then...

    # move to the snakemake's working directory
    cd cola-seq_bp_calling/code
    # Create environment for the snakemake
    mamba env create -f envs/cola-seq_bp_calling.yaml
    # And activate the enviroment
    conda activate cola-seq_bp_calling

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config/config.yaml`. Use/modify the config yaml files in the `snakemake_profiles/slurm/` profile to run on UChicago RCC Midway with slurm scheduler. Define samples tsv file referenced in `config/config.yaml` (eg `.test/samples.tsv` by default). See comments in `config/config.yaml` for more description of configuration. See below for description of the sample tsv file:

#### samples tsv file:

Comment lines (prefixed iwth '#') are ignored.  Each line must contain unique sample name.  If filepath is provided in Bam_file is and the file exists, the snakemake pipeline will use that. Otherwise, will align from paired gzipped fastq using rules defined in Snakefile (see 'include' statements). Fastq and bam must be paired end. If a sample has multiple fastq (eg a single library sequenced on two sequencing lanes creates two fastq), separate multiple fastq files by commas (no spaces) in the appropriate cell. Samples with multiple fastq will get combined before alignment to output a single bam.  Samples with the same BP_group value will get merged at the level of bam files before performing BP_calling. This test data demonstrates the pipeline. The first sample references a bam. The rest of the samples reference gzipped fastq files. The last two samples have two paired fastq files per sample.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via the included slurm snakemake profile.

    snakemake --profile snakemake_profiles/slurm

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
