rule STAR_make_index:
    """
    did not work on bigmem2. Never figured out why (the log file didn't
    indicate anything). Ran on login node with success.
    """
    input:
        fasta=config["fasta"],
        gtf=config["gtf"],
    output:
        index="ReferenceGenome/STARIndex/chrLength.txt",
    log:
        "logs/STAR_make_index.log",
    params:
        genomeDir="ReferenceGenome/STARIndex/",
    threads: 4
    resources:
        mem="42G",
        partition="bigmem2",
        ntasks=5,
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN {threads} --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """


rule CopyAndMergeFastq:
    input:
        R1=GetFastqList_For_Sample("Fastq_R1_list"),
        R2=GetFastqList_For_Sample("Fastq_R2_list"),
    output:
        R1="Fastq/{sample}/R1.fastq.gz",
        R2="Fastq/{sample}/R2.fastq.gz",
    log:
        "logs/CopyAndMergeFastq/{sample}.log",
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log}
        cat {input.R2} > {output.R2} 2>> {log}
        """


rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1="Fastq/{sample}/R1.fastq.gz",
        R2="Fastq/{sample}/R2.fastq.gz",
    output:
        R1="FastqFastp/{sample}/R1.fastq.gz",
        R2="FastqFastp/{sample}/R2.fastq.gz",
        html="FastqFastp/{sample}/fastp.html",
        json="FastqFastp/{sample}/fastp.json",
    params:
        Get_fastq_params,
    resources:
        mem_mb=8000,
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} {params} --html {output.html} --json {output.json} &> {log}
        """


rule STAR_Align:
    input:
        index="ReferenceGenome/STARIndex/chrLength.txt",
        R1="FastqFastp/{sample}/R1.fastq.gz",
        R2="FastqFastp/{sample}/R2.fastq.gz",
    output:
        bam="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
        align_log="Alignments/STAR_Align/{sample}/Log.final.out",
    threads: 8
    log:
        "logs/STAR_Align/{sample}.log",
    params:
        readMapNumber=-1,
        ENCODE_params="--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
    resources:
        cpus_per_node=9,
        mem=58000,
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix Alignments/STAR_Align/{wildcards.sample}/ --genomeDir ReferenceGenome/STARIndex/ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 --limitBAMsortRAM 16000000000 {params.ENCODE_params} --outSAMstrandField intronMotif  &> {log}
        samtools index {output.bam}
        """


rule DedupStarAlignments:
    """
    If sample is from test config, it already has been deduped, so just make a
    mock bam file with cp
    """
    input:
        bam="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        bam=temp("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.dedup.bam"),
    resources:
        mem_mb=12000,
    log:
        "logs/DedupStarAlignments/{sample}.log",
    conda:
        "../envs/umi_tools.yml"
    params:
        IsSampleFromTestConfig_smk_input_function,
    shell:
        """
        if {params}
        then
            cp {input.bam} {output.bam}
        else
            umi_tools dedup --paired -I {input.bam} -S {output.bam} -L {log} --umi-separator=:
        fi
        """
