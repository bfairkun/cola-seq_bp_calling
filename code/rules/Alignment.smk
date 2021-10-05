rule STAR_make_index:
    """
    did not work on bigmem2. Never figured out why (the log file didn't
    indicate anything). Ran on login node with success.
    """
    input:
        fasta = config["fasta"]
        gtf = config["gtf"],
    output:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
    log:
        "logs/STAR_make_index.log"
    params:
        genomeDir = "ReferenceGenome/STARIndex/"
    threads: 4
    resources:
        mem = "42G",
        partition = "bigmem2",
        ntasks = 5
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN {threads} --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule CopyAndMergeFastq:
    input:
    output:
        output
    log:
        "logs/CopyAndMergeFastq.log"
    shell:
        """
        shell
        """

rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "Fastq/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz"
    output:
        R1 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz",
        html = "FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.html",
        json = "FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.json"
    params:
        umi = GetFastpParamsUmi2,
        I = "-I",
        O = "-O"
    resources:
        mem_mb = 8000
    log:
        "logs/fastp/{Phenotype}.{IndID}.{Rep}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} {params.I} {input.R2} -o {output.R1} {params.O} {output.R2} --html
        """

rule STAR_Align_WASP:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz",
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf"
    output:
        bam = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam",
        align_log = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align_WASP/{Phenotype}/{IndID}.{Rep}.log"
    params:
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        WASP_params = "--waspOutputMode SAMtag --outSAMattributes NH HI AS nM XS vW --varVCFfile"
    resources:
        cpus_per_node = 9,
        mem = 58000,
    wildcard_constraints:
        Phenotype = "Expression.Splicing|chRNA.Expression.Splicing"
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix Alignments/STAR_Align/{wildcards.Phenotype}/{wildcards.IndID}/{wildcards.Rep}/ --genomeDir ReferenceGenome/STARIndex/ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 {params.WASP_params} {input.vcf} --limitBAMsortRAM 16000000000 {params.ENCODE_params} --outSAMstrandField intronMotif  &> {log}
        """

rule DedupStarAlignments:
    input:
        bam = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam.ba
    output:
        bam = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.dedup.
    resources:
        mem_mb = 12000
    log:
        "logs/DedupStarAlignments/{Phenotype}/{IndID}.{Rep}.log"
    conda:
        "../envs/umi_tools.yml"
    params:
        GetUmitoolsDedupParams
    shell:
        """
        umi_tools dedup --paired -I {input.bam} -S {output.bam} -L {log} --umi-separator=: {par
        """
