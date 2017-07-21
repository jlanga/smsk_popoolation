rule mpileup_population_chromosome:
    input:
        cram = MAP_FILT + "{population}/{chromosome}.cram",
        fa  = RAW + "genome.fa",
        fai = RAW + "genome.fa.fai"
    output:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz"
    threads:
        1
    log:
        MPILEUP_RAW + "{population}/{chromosome}.log"
    benchmark:
        MPILEUP_RAW + "{population}/{chromosome}.json"
    shell:
        "(samtools mpileup "
            "-B "
            "-Q 0 "
            "-f {input.fa} "
            "{input.cram} "
            "| pigz "
                "--best ) "
        "> {output.mpileup_gz} "
        "2> {log}"



rule mpileup_filter_population_chromosome_gtf:
    """
    Get a GTF with the indels present
    """
    input:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz"
    output:
        gtf = temp(
            MPILEUP_FILT + "{population}/{chromosome}.gtf"
        )
    threads:
        1
    params:
        indel_window = config["popoolation_params"]["find_indels"]["indel_window"],
        min_count    = config["popoolation_params"]["find_indels"]["min_count"]
    log:
        MPILEUP_FILT + "{population}/{chromosome}.gtf.log"
    benchmark:
        MPILEUP_FILT + "{population}/{chromosome}.gtf.json"
    shell:
        "perl src/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl "
            "--input <( pigz --decompress --stdout {input.mpileup_gz} ) "
            "--output {output.gtf} "
            "--indel-window {params.indel_window} "
            "--min-count {params.min_count} "
        "2>{log} 1>&2"



rule mpileup_filter_population_chromosome_mpileup:
    """
    Filter indels from an mpileup given a GTF files with their coordinates
    """
    input:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz",
        gtf = MPILEUP_FILT + "{population}/{chromosome}.gtf"
    output:
        mpileup_gz = temp(
            MPILEUP_FILT + "{population}/{chromosome}.mpileup.gz"
        )
    params:
        mpileup = MPILEUP_FILT + "{population}/{chromosome}.mpileup"
    threads:
        1
    log:
        MPILEUP_FILT + "{population}/{chromosome}.mpileup.log"
    benchmark:
        MPILEUP_FILT + "{population}/{chromosome}.mpileup.json"
    shell:
        "perl src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl "
            "--input <(pigz --decompress --stdout {input.mpileup_gz} ) "
            "--gtf {input.gtf} "
            "--output {params.mpileup} "
        "2> {log} 1>&2 ; "
        "pigz --best {params.mpileup} 2>> {log}"



rule mpileup_subsample_population_chromosome:
    input:
        mpileup_gz = MPILEUP_FILT + "{population}/{chromosome}.mpileup.gz"
    output:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    params:
        mpileup        = MPILEUP_SUB + "{population}/{chromosome}.mpileup",
        minqual        = config["popoolation_params"]["subsample"]["minqual"],
        method         = config["popoolation_params"]["subsample"]["method"],
        maxcoverage    = config["popoolation_params"]["subsample"]["maxcoverage"],
        targetcoverage = config["popoolation_params"]["subsample"]["targetcoverage"]
    threads:
        1
    log:
        MPILEUP_SUB + "{population}/{chromosome}.log"
    benchmark:
        MPILEUP_SUB + "{population}/{chromosome}.json"
    shell:
        "perl src/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl "
            "--min-qual {params.minqual} "
            "--method {params.method} "
            "--max-coverage {params.maxcoverage} "
            "--fastq-type sanger "
            "--target-coverage {params.targetcoverage} "
            "--input <(pigz --decompress --stdout {input.mpileup_gz}) "
            "--output {params.mpileup} "
        "2> {log} 1>&2 ; "
        "pigz --best {params.mpileup} 2>> {log}"
