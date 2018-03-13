rule map_bwa_index:
    """Index with bwa"""
    input:
        fa = RAW + "genome.fa"
    output:
        mock  = touch(
            MAP_RAW + "genome"
        ),
        buckets = expand(
            MAP_RAW + "genome.{suffix}",
            suffix = "amb ann bwt pac sa".split()
        )
    log: MAP_RAW + "genome.bwa_index.log"
    benchmark: MAP_RAW + "genome.bwa_index.json"
    conda: "map.yml"
    shell:
        "bwa index -p {output.mock} {input.fa} > {log} 2>&1"



rule map_bwa_map:
    """Map population with bowtie2, sort with samtools, compress to cram"""
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
        bowtie2_params = config["bwa_params"],
        sq_id = "{population}_{library}",
        sq_library = "LB:truseq_{library}",
        sq_platform = "PL:Illumina",
        sq_sample = "SM:{population}",
    threads: 24
    log: MAP_RAW + "{population}/{library}.bwa_mem.log"
    benchmark: MAP_RAW + "{population}/{library}.bwa_mem.json"
    conda: "map.yml"
    shell:
        "(bwa mem "
            "-M "
            "-R '@RG\tID:{params.sq_id}\t{params.sq_library}\t{params.sq_platform}\t{params.sq_sample}' "
            "-t {threads} "
            "{params.bowtie2_params} "
            "{input.index} "
            "{input.forward} "
            "{input.reverse} "
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
