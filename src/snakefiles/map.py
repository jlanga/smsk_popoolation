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
    input:
        forward = QC + "{population}_1.fq.gz",
        reverse = QC + "{population}_2.fq.gz",
        index = MAP_RAW + "genome"
    output:
        bam = protected(
            MAP_RAW + "{population}.bam"
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
        "fifo_aux=$(mktemp --dry-run) ; "
        "mkfifo $fifo_aux ; "
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
            "> $fifo_aux "
            "| picard SortSam "
                "INPUT=$fifo_aux "
                "OUTPUT={output.bam} "
                "COMPRESSION_LEVEL=9 "
                "VALIDATION_STRINGENCY=SILENT "
                "SORT_ORDER=coordinate ) "
        "2> {log}"



rule map_population_bai:
    input:
        bam = MAP_RAW + "{population}.bam"
    output:
        bai = MAP_RAW + "{population}.bam.bai"
    threads:
        1
    log:
        MAP_RAW + "{population}.bai.log"
    benchmark:
        MAP_RAW + "{population}.bai.json"
    shell:
        "samtools index {input.bam} 2> {log}"



rule map_split_population_chromosome_split:
    """
    We use uncompressed bam to accelerate the output. The result of this rule is temporary
    """
    input:
        bam = MAP_RAW + "{population}.bam",
        bai = MAP_RAW + "{population}.bam.bai"
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
            "{input.bam} "
            "{params.chromosome} "
        "> {output.bam} "
        "2> {log}"



rule map_filter_population_chromosome:
    input:
        bam = MAP_SPLIT + "{population}/{chromosome}.bam",
    output:
        bam = protected(
            MAP_FILT + "{population}/{chromosome}.bam"
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
        "fifo_aux1=$(mktemp --dry-run) ; mkfifo $fifo_aux1 ; "
        "fifo_aux2=$(mktemp --dry-run) ; mkfifo $fifo_aux2 ; "
        "( picard -Xmx4g MarkDuplicates "
            "INPUT={input.bam} "
            "OUTPUT=$fifo_aux1 "
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
                "$fifo_aux1 "
                "> $fifo_aux2 "
            "| picard SortSam "
                "INPUT=$fifo_aux2 "
                "OUTPUT={output.bam} "
                "VALIDATION_STRINGENCY=SILENT "
                "SORT_ORDER=coordinate "
                "COMPRESSION_LEVEL=9 ) "
        "2> {log} ; "
        "rm $fifo_aux1 $fifo_aux2"    
