rule map_bowtie2_index:
    """Index with bowtie2"""
    input:
        fa = RAW + "genome.fa"
    output:
        mock  = touch(
            MAP_RAW + "genome"
        ),
        buckets = expand(
            MAP_RAW + "genome.{suffix}.bt2",
            suffix = "1 2 3 4 rev.1 rev.2".split()
        )
    log: MAP_RAW + "genome.bowtie2_index.log"
    benchmark: MAP_RAW + "genome.bowtie2_index.json"
    conda: "map.yml"
    shell:
        "bowtie2-build {input.fa} {output.mock} > {log} 2>&1"



rule map_bowtie2_map:
    """Map population with bowtie2, sort with picard, compress to cram with samtools"""
    input:
        forward = QC + "{population}/{library}_1.fq.gz",
        reverse = QC + "{population}/{library}_2.fq.gz",
        unp_forward = QC + "{population}/{library}_3.fq.gz",
        unp_reverse = QC + "{population}/{library}_4.fq.gz",
        index = MAP_RAW + "genome",
        reference = RAW + "genome.fa"
    output:
        cram = protected(
            MAP_RAW + "{population}/{library}.cram"
        )
    params:
        bowtie2_params = config["bowtie2_params"],
        sq_id = "{population}_{library}",
        sq_library = "LB:truseq_{library}",
        sq_platform = "PL:Illumina",
        sq_sample = "SM:{population}",
    threads: 24
    log: MAP_RAW + "{population}/{library}.bowtie2.log"
    benchmark: MAP_RAW + "{population}/{library}.bowtie2.json"
    conda: "map.yml"
    shell:
        "(bowtie2 "
            "--rg-id {params.sq_id} "
            "--rg {params.sq_library} "
            "--rg {params.sq_platform} "
            "--rg {params.sq_sample} "
            "--threads {threads} "
            "{params.bowtie2_params} "
            "-x {input.index} "
            "-1 {input.forward} "
            "-2 {input.reverse} "
            "-U {input.unp_forward},{input.unp_reverse} "
        "| samtools sort "
            "-l 9 "
            "-o {output.cram} "
            "--reference {input.reference} "
            "--output-fmt CRAM "
            "-@ {threads} "
            "/dev/stdin "
        ") 2> {log}"



rule map_split:  # USE BAM bc it markduplicates needs a file
    """Extract chromosome in cram

    We use uncompressed bam to accelerate the output. The result of this rule is
    temporary.

    Note: the following step is picard MarkDuplicates, and needs a proper file
    since it makes two passes.
    """
    input:
        cram = MAP_RAW + "{population}/{library}.cram",
        crai = MAP_RAW + "{population}/{library}.cram.crai",
        reference = RAW + "genome.fa"
    output:
        bam = temp(
            MAP_SPLIT + "{population}/{library}/{chromosome}.bam"
        )
    threads: 1
    params:
        chromosome = "{chromosome}"
    log: MAP_SPLIT + "{population}/{library}/{chromosome}.log"
    benchmark: MAP_SPLIT + "{population}/{library}/{chromosome}.json"
    conda: "map.yml"
    shell:
        "samtools view "
            "-u "
            "-T {input.reference} "
            "-o {output.bam} "
            "-@ {threads} "
            "{input.cram} "
            "{params.chromosome} "
        "2> {log}"



rule map_filter:  # TODO: java memory, uncompressed bam
    """Remove duplicates from CRAM and filter out sequences.

    samtools view | MarkDuplicates | samtools view -f -F | SortSam
    samtools view

    Pairs with something unpaired will disappear.
    """
    input:
        bam = MAP_SPLIT + "{population}/{library}/{chromosome}.bam",
        reference = RAW + "genome.fa"
    output:
        cram = protected(
            MAP_FILT + "{population}/{library}/{chromosome}.cram"
        ),
        dupstats = MAP_FILT + "{population}/{library}/{chromosome}.dupstats"
    params:
        memory= config["picard_markduplicates_params"]["memory"]
    log: MAP_FILT + "{population}/{library}/{chromosome}.log"
    benchmark: MAP_FILT + "{population}/{library}/{chromosome}.json"
    threads: 1
    conda: "map.yml"
    shell:
        "(picard -Xmx{params.memory} MarkDuplicates "
            "INPUT={input.bam} "
            "OUTPUT=/dev/stdout "
            "METRICS_FILE={output.dupstats} "
            "ASSUME_SORT_ORDER=coordinate "
            "VALIDATION_STRINGENCY=SILENT "
            "COMPRESSION_LEVEL=0 "
            "REMOVE_DUPLICATES=true "
            "QUIET=false "
        "| samtools view "
            "-q 20 "
            "-f 0x0002 "  # read mapped in proper pair. Leave only
            "-F 0x0004 "  # read unmapped. Throw away
            "-F 0x0008 "  # mate unmapped. Throw away
            "-u "
            "- "
        "| samtools sort "
            "-l 9 "
            "-o {output.cram} "
            "--reference {input.reference} "
            "--output-fmt CRAM "
            "-@ {threads} "
            "/dev/stdin "
        ") 2> {log}"
