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


rule preprocess__map__split__:
    """Extract chromosome in cram

    We use uncompressed bam to accelerate the output. The result of this rule
    is temporary.

    Note: the following step is picard MarkDuplicates, and needs a proper file
    since it makes two passes. Output is a bam because MarkDuplicates needs
    one.
    """
    input:
        cram=PRE_MAP / "{population}.{library}.cram",
        crai=PRE_MAP / "{population}.{library}.cram.crai",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        bam=temp(PRE_SPLIT / "{population}.{library}.{chromosome}.bam"),
    params:
        chromosome="{chromosome}",
    log:
        PRE_SPLIT / "{population}.{library}.{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        samtools view \
            -u \
            -T {input.reference} \
            -o {output.bam} \
            {input.cram} \
            {params.chromosome} \
        2> {log}
        """


rule preprocess__map__filter__:  # TODO: java memory, uncompressed bam
    """Remove duplicates from CRAM and filter out sequences.

    samtools view | MarkDuplicates | samtools view -f -F | SortSam
    samtools view

    Pairs with something unpaired will disappear.

    -f 0x0002  # read mapped in proper pair. Leave only
    -F 0x0004  # read unmapped. Throw away
    -F 0x0008  # mate unmapped. Throw away
    """
    input:
        bam=PRE_SPLIT / "{population}.{library}.{chromosome}.bam",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        cram=PRE_FILT / "{population}.{library}.{chromosome}.cram",
        dupstats=PRE_FILT / "{population}.{library}.{chromosome}.dupstats",
    log:
        PRE_FILT / "{population}.{library}.{chromosome}.log",
    resources:
        memory_gb=params["picard_markduplicates"]["memory_gb"],
    conda:
        "__environment__.yml"
    shell:
        """
        ( picard MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT /dev/stdout \
            --METRICS_FILE {output.dupstats} \
            --ASSUME_SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY SILENT \
            --COMPRESSION_LEVEL 0 \
            --REMOVE_DUPLICATES true \
            --QUIET true \
            --VERBOSITY ERROR \
        | samtools view \
            -q 20 \
            -f 2 \
            -F 4 \
            -F 8 \
            -@ {threads} \
            --output-fmt cram \
            --reference {input.reference} \
            -o {output.cram} \
            - \
        ) 2> {log}
        """


rule preprocess__map:
    input:
        [
            PRE_FILT / f"{population}.{library}.{chromosome}.cram"
            for population, library in POPULATION_LIBRARY
            for chromosome in CHROMOSOMES
        ],
