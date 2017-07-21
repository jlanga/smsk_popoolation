rule map_genome_bowtie2_index:
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
    threads:
        1
    log:
        MAP_RAW + "genome.bowtie2_index.log"
    benchmark:
        MAP_RAW + "genome.bowtie2_index.json"
    shell:
        "bowtie2-build "
            "{input.fa} "
            "{output.mock} "
        "> {log} 2>&1"



rule map_population_bowtie2:
    """
    Map population with bowtie2, sort with picard, compress to cram with
    samtools.
    """
    input:
        forward = QC + "{population}_1.fq.gz",
        reverse = QC + "{population}_2.fq.gz",
        index = MAP_RAW + "genome",
        reference = RAW + "genome.fa"
    output:
        cram = protected(
            MAP_RAW + "{population}.cram"
        )
    params:
        bowtie2_params = config["bowtie2_params"],
        sq_id = "{population}",
        sq_library = "LB:truseq_{population}",
        sq_platform = "PL:Illumina",
        sq_sample = "SM:{population}",
    threads:
        24
    log:
        MAP_RAW + "{population}.bowtie2.log"
    benchmark:
        MAP_RAW + "{population}.bowtie2.json"
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
        "| picard SortSam "
            "INPUT=/dev/stdin "
            "OUTPUT=/dev/stdout "
            "COMPRESSION_LEVEL=0 "
            "VALIDATION_STRINGENCY=SILENT "
            "SORT_ORDER=coordinate "
        "| samtools view "
            "-@ {threads} "
            "-T {input.reference} "
            "-C "
            "-o {output.cram} "
            "/dev/stdin"
        ") 2> {log}"



rule map_split_population_chromosome_split:  # USE BAM bc it markduplicates needs a file
    """
    We use uncompressed bam to accelerate the output. The result of this rule is
    temporary.

    Note: the following step is picard MarkDuplicates, and needs a proper file
    since it makes two passes over it
    """
    input:
        cram = MAP_RAW + "{population}.cram",
        crai = MAP_RAW + "{population}.cram.crai",
        reference = RAW + "genome.fa"
    output:
        bam = temp(
            MAP_SPLIT + "{population}/{chromosome}.bam"
        )
    params:
        chromosome = "{chromosome}"
    threads:
        1
    log:
        MAP_SPLIT + "{population}/{chromosome}.log"
    benchmark:
        MAP_SPLIT + "{population}/{chromosome}.json"
    shell:
        "samtools view "
            "-u "
            "-T {input.reference} "
            "-o {output.bam} "
            "{input.cram} "
            "{params.chromosome} "
        "2> {log}"



rule map_filter_population_chromosome:  # TODO: java memory
    """
    Remove duplicates from CRAM and filter out sequences.

    samtools view | MarkDuplicates | samtools view -f -F | SortSam | \
    samtools view
    """
    input:
        bam = MAP_SPLIT + "{population}/{chromosome}.bam",
        # crai = MAP_RAW + "{population}.cram.crai",
        reference = RAW + "genome.fa"
    output:
        cram = protected(
            MAP_FILT + "{population}/{chromosome}.cram"
        ),
        dupstats = MAP_FILT + "{population}/{chromosome}.dupstats"
    params:
        chromosome = "{chromosome}"
    log:
        MAP_FILT + "{population}/{chromosome}.log"
    benchmark:
        MAP_FILT + "{population}/{chromosome}.json"
    threads:
        4 # Too much ram usage
    shell:
        "(picard -Xmx4g MarkDuplicates "
            "INPUT={input.bam} "
            "OUTPUT=/dev/stdout "
            "METRICS_FILE={output.dupstats} "
            "ASSUME_SORTED=true "
            "VALIDATION_STRINGENCY=SILENT "
            "COMPRESSION_LEVEL=0 "
            "REMOVE_DUPLICATES=true "
            "QUIET=false "
        "| samtools view "
            "-q 20 "
            "-f 0x0002 "
            "-F 0x0004 "
            "-F 0x0008 "
            "-u "
            "- "
        "| picard SortSam "
            "INPUT=/dev/stdin "
            "OUTPUT=/dev/stdout "
            "VALIDATION_STRINGENCY=SILENT "
            "SORT_ORDER=coordinate "
            "COMPRESSION_LEVEL=0 "
        "| samtools view "
            "-C "
            "-T {input.reference} "
            "-o {output.cram} "
            "/dev/stdin "
        ") 2> {log}"
