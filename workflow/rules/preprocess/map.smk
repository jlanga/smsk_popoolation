rule preprocess__map__bwamem2__:
    """Map population with bowtie2, sort with samtools, compress to cram"""
    input:
        forward_=READS / "{population}.{library}_1.fq.gz",
        reverse_=READS / "{population}.{library}_2.fq.gz",
        index=PRE_INDEX / f"{REFERENCE_NAME}",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        cram=PRE_MAP / "{population}.{library}.cram",
    params:
        rg_tag=compose_rg_tag,
    threads: 24
    log:
        PRE_MAP / "{population}.{library}.bwa_mem.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( bwa-mem2 mem \
            -M \
            -R '{params.rg_tag}' \
            -t {threads} \
            {input.index} \
            {input.forward_} \
            {input.reverse_} \
        | samtools sort \
            -o {output.cram} \
            --reference {input.reference} \
            --output-fmt CRAM \
            -@ {threads} \
        ) 2> {log}
        """


rule preprocess__map__rmdup__:
    """Run GATK MarkDuplicates to find AND remove duplicates"""
    input:
        cram=PRE_MAP / "{population}.{library}.cram",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        bam=temp(PRE_RMDUP / "{population}.{library}.bam"),
        dupstats=PRE_RMDUP / "{population}.{library}.dupstats",
    log:
        PRE_RMDUP / "{population}.{library}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        gatk MarkDuplicates \
            --INPUT {input.cram} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.dupstats} \
            --REFERENCE_SEQUENCE {input.reference} \
            --ASSUME_SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY SILENT \
            --COMPRESSION_LEVEL 0 \
            --REMOVE_DUPLICATES true \
            --QUIET true \
            --VERBOSITY ERROR \
        2> {log}
        """


rule preprocess__map__filter__:
    """Remove unpaired and unmapped sequences

    -f 0x0002  : read mapped and in proper pair. Keep
    -F 0x0004  : read unmapped. Throw away
    -F 0x0008  : mate unmapped. Throw away
    """
    input:
        bam=PRE_RMDUP / "{population}.{library}.bam",
    output:
        bam=temp(PRE_FILT / "{population}.{library}.bam"),
    log:
        PRE_FILT / "{population}.{library}.log",
    resources:
        memory_gb=params["picard_markduplicates"]["memory_gb"],
    conda:
        "__environment__.yml"
    shell:
        """
        samtools view \
            --min-MQ 20 \
            --require-flags 2 \
            --exclude-flags 4 \
            --exclude-flags 8 \
            --uncompressed \
            --output {output.bam} \
            {input.bam} \
        2> {log}
        """


rule preprocess__map__split__:
    """Extract only one chromosome and store it in a cram file
    """
    input:
        bam=PRE_FILT / "{population}.{library}.bam",
        bai=PRE_FILT / "{population}.{library}.bam.bai",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        cram=PRE_SPLIT / "{population}.{library}" / "{chromosome}.cram",
    params:
        chromosome=lambda w: w.chromosome,
    log:
        PRE_SPLIT / "{population}.{library}" / "{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        samtools view \
            --uncompressed \
            --reference {input.reference} \
            --output {output.cram} \
            {input.bam} \
            {params.chromosome} \
        2> {log}
        """


rule preprocess__map:
    input:
        [
            PRE_SPLIT / f"{population}.{library}" / f"{chromosome}.cram"
            for population, library in POPULATION_LIBRARY
            for chromosome in CHROMOSOMES
        ],
