rule preprocess__map__index:
    """Index with bwa"""
    input:
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        mock=touch(MAP_INDEX / f"{REFERENCE_NAME}"),
        buckets=[
            MAP_INDEX / f"{reference_name}.{suffix}"
            for reference_name in [REFERENCE_NAME]
            for suffix in "amb ann bwt pac sa".split()
        ],
    log:
        MAP_INDEX / "bwa_index.log",
    conda:
        "__environment__.yml"
    shell:
        "bwa index -p {output.mock} {input.fa} > {log} 2>&1"


rule preprocess__map__bwamem__:
    """Map population with bowtie2, sort with samtools, compress to cram"""
    input:
        forward_=READS / "{population}.{library}_1.fq.gz",
        reverse_=READS / "{population}.{library}_2.fq.gz",
        index=MAP_INDEX / f"{REFERENCE_NAME}",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        cram=MAP_RAW / "{population}.{library}.cram",
    params:
        extra=params["bwa"]["extra"],
        rg_tag=compose_rg_tag,
    threads: 24
    log:
        MAP_RAW / "{population}.{library}.bwa_mem.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( bwa mem \
            -M \
            -R '{params.rg_tag}' \
            -t {threads} \
            {params.extra} \
            {input.index} \
            {input.forward_} \
            {input.reverse_} \
        | samtools sort \
            -l 9 \
            -o {output.cram} \
            --reference {input.reference} \
            --output-fmt CRAM \
            -M \
            -@ {threads} \
            `#/dev/stdin` \
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
        cram=MAP_RAW / "{population}.{library}.cram",
        crai=MAP_RAW / "{population}.{library}.cram.crai",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        bam=temp(MAP_SPLIT / "{population}.{library}.{chromosome}.bam"),
    params:
        chromosome="{chromosome}",
    log:
        MAP_SPLIT / "{population}.{library}.{chromosome}.log",
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
    """
    input:
        bam=MAP_SPLIT / "{population}.{library}.{chromosome}.bam",
        reference=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        cram=MAP_FILT / "{population}.{library}.{chromosome}.cram",
        dupstats=MAP_FILT / "{population}.{library}.{chromosome}.dupstats",
    log:
        MAP_FILT / "{population}.{library}.{chromosome}.log",
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
            -f 0x0002  `# read mapped in proper pair. Leave only` \
            -F 0x0004  `# read unmapped. Throw away` \
            -F 0x0008  `# mate unmapped. Throw away` \
            -u \
            - \
        | samtools sort \
            -l 9 \
            -@ {threads} \
            -o {output.cram} \
            --reference {input.reference} \
            --output-fmt CRAM \
            /dev/stdin \
        ) 2> {log}
        """


rule preprocess__map:
    input:
        [
            MAP_FILT / f"{population}.{library}.{chromosome}.cram"
            for population, library in POPULATION_LIBRARY
            for chromosome in CHROMOSOMES
        ],
