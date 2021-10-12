# test data description

test data to run snakemake pipeline on. All of the data is in some form data from Yi's CoLa-seq libraries around the URB1 gene (which is pictured in Fig1D). See commented desriptions of files in tree below to understand what the files are.

```
.test/test_data/
├── README.md
    # All the files starting with Split* are a subset of reads from
    # testdata.bazam.R[12].fastq.gz.
    # Therefore (see below) they contain all the reads in the specified region
    # and their mate pairs. The subsets are based on separate libraries that Yi
    # prepared.
├── Split_camptothecin.YZ132.R1.fastq.gz
├── Split_camptothecin.YZ132.R2.fastq.gz
├── Split_DMSO_ctrl_for_PlaB.R1.fastq.gz
├── Split_DMSO_ctrl_for_PlaB.R2.fastq.gz
├── Split_DMSO.YZ131.R1.fastq.gz
├── Split_DMSO.YZ131.R2.fastq.gz
├── Split_flavopiridol.YZ133.R1.fastq.gz
├── Split_flavopiridol.YZ133.R2.fastq.gz
├── Split_PlaB.R1.fastq.gz
├── Split_PlaB.R2.fastq.gz
├── Split_YZ111.R1.fastq.gz
├── Split_YZ111.R2.fastq.gz
├── Split_YZ112.R1.fastq.gz
├── Split_YZ112.R2.fastq.gz
├── Split_YZ151_1.R1.fastq.gz
├── Split_YZ151_1.R2.fastq.gz
├── Split_YZ151_2.R1.fastq.gz
├── Split_YZ151_2.R2.fastq.gz
├── Split_YZ152_1.R1.fastq.gz
├── Split_YZ152_1.R2.fastq.gz
├── Split_YZ152_2.R1.fastq.gz
├── Split_YZ152_2.R2.fastq.gz
├── testdata.bam
    # the bamfile output of:
    # samtools view -bh bams/By/Chrom/Merged.chr21.bam chr21:32,309,018-32,395,012
    # where bams/By/Chrom/Merged.chr21.bam is a merged bam of all Yi's samples
    # therefore, reads not in the specified region are not in the bam, even if
    # their mate pair is
├── testdata.bam.bai
├── testdata.bazam.R1.fastq.gz
└── testdata.bazam.R2.fastq.gz
    # The R1 and R2 fastq output of:
    # bazam -bam Merged.chr21.bam -L chr21:32,309,018-32,395,012 -r2 ~/cola-seq_bp_calling/code/.test/test_data/testdata.bazam.R1.fastq -r1 ~/cola-seq_bp_calling/code/.test/test_data/testdata.bazam.R2.fastq
    # Therefore, these contain all the reads in the specied region and their mate pairs
    # Note that I swapped r1 and r2. This is because there is an issue with
    # bazam (https://github.com/ssadedin/bazam/issues/31) wherein read pairs
    # where R1 is on the minus strand (such as the minus stranded gene URB1 that
    # the test data is based on) get swapped by bazam

0 directories, 27 files
```
