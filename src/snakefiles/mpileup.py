rule mpileup_population_chromosome:
    """
    Compute the mpileup and compress it
    """
    input:
        cram = MAP_FILT + "{population}/{chromosome}.cram",
        fa  = RAW + "genome.fa",
        fai = RAW + "genome.fa.fai"
    output:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz"
    log: MPILEUP_RAW + "{population}/{chromosome}.log"
    benchmark: MPILEUP_RAW + "{population}/{chromosome}.json"
    threads: 2  # mpileup and gzip work at the same pace
    shell:
        "(samtools mpileup "
            "--no-BAQ "
            "--min-BQ 0 "
            "--fasta-ref {input.fa} "
            "{input.cram} "
        "| gzip "
            "--best "
        "> {output.mpileup_gz} "
        ") 2> {log}"



rule mpileup_filter_population_chromosome_gtf:
    """
    Get a GTF with the indels present.
    """
    input:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz"
    output:
        gtf = temp(
            MPILEUP_FILT + "{population}/{chromosome}.gtf"
        )
    params:
        indel_window = config["popoolation_params"]["find_indels"]["indel_window"],
        min_count    = config["popoolation_params"]["find_indels"]["min_count"]
    log: MPILEUP_FILT + "{population}/{chromosome}.gtf.log"
    benchmark: MPILEUP_FILT + "{population}/{chromosome}.gtf.json"
    shell:
        "perl src/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl "
            "--input <(gzip --decompress --stdout {input.mpileup_gz}) "
            "--output {output.gtf} "
            "--indel-window {params.indel_window} "
            "--min-count {params.min_count} "
        "2>{log} 1>&2"



rule mpileup_filter_population_chromosome_mpileup:
    """
    Filter indels from an mpileup given a GTF files with their coordinates.
    Compress results.
    """
    input:
        mpileup_gz = MPILEUP_RAW + "{population}/{chromosome}.mpileup.gz",
        gtf = MPILEUP_FILT + "{population}/{chromosome}.gtf"
    output:
        mpileup_fifo = temp(
            MPILEUP_FILT + "{population}/{chromosome}.mpileup"
        ),
        mpileup_gz = temp(
            MPILEUP_FILT + "{population}/{chromosome}.mpileup.gz"
        )
    log: MPILEUP_FILT + "{population}/{chromosome}.mpileup.log"
    benchmark: MPILEUP_FILT + "{population}/{chromosome}.mpileup.json"
    shell:
        "mkfifo {output.mpileup_fifo}; "
        "(cat {output.mpileup_fifo} | gzip --fast > {output.mpileup_gz} &);"
        "perl src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl "
            "--input <(gzip --decompress --stdout {input.mpileup_gz}) "
            "--gtf {input.gtf} "
            "--output {output.mpileup_fifo} "
        "2> {log} 1>&2"



rule mpileup_subsample_population_chromosome:
    """
    Perform the subsampling step. Compress results as they are generated through
    a FIFO.
    """
    input:
        mpileup = MPILEUP_FILT + "{population}/{chromosome}.mpileup.gz"
    output:
        mpileup_fifo = temp(
            MPILEUP_SUB + "{population}/{chromosome}.mpileup"
        ),
        mpileup_gz = protected(
            MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
        )
    params:
        minqual        = config["popoolation_params"]["subsample"]["minqual"],
        method         = config["popoolation_params"]["subsample"]["method"],
        maxcoverage    = config["popoolation_params"]["subsample"]["maxcoverage"],
        targetcoverage = config["popoolation_params"]["subsample"]["targetcoverage"]
    log: MPILEUP_SUB + "{population}/{chromosome}.log"
    benchmark: MPILEUP_SUB + "{population}/{chromosome}.json"
    shell:
        "mkfifo {output.mpileup_fifo}; "
        "(cat {output.mpileup_fifo} | gzip --best > {output.mpileup_gz} &); "
        "perl src/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl "
            "--min-qual {params.minqual} "
            "--method {params.method} "
            "--max-coverage {params.maxcoverage} "
            "--fastq-type sanger "
            "--target-coverage {params.targetcoverage} "
            "--input <(gzip -dc {input.mpileup}) "
            "--output {output.mpileup_fifo} "
        "2> {log} 1>&2 ; "
