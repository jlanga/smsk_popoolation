rule sync_identify_indels:
    """Identify indels like in mpileup_popoolation_identify_indels, but all together"""
    input:
        mpileups = expand(
            MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
            population=POPULATIONS,
            chromosome="{chromosome}"
        )
    output:
        gtf = temp(SYNC_FILT + "{chromosome}.gtf")
    params:
        indel_window = config["popoolation2_params"]["find_indels"]["indel_window"],
        min_count    = config["popoolation2_params"]["find_indels"]["min_count"],
        mpileups_comma = "{" + ",".join(
            expand(
                MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
                population=POPULATIONS,
                chromosome="{chromosome}"
            )
        ) + "}",
    threads: 1
    log: SYNC_FILT + "{chromosome}.identify_indels.log"
    benchmark: SYNC_FILT + "{chromosome}.identify_indels.bmk"
    conda: "sync.yml"
    shell:
        "(eval "
            "paste <(gzip -dc {input.mpileups[0]} | cut -f 1-3) "
            "'<(gzip -dc '{params.mpileups_comma}' | cut -f 4-6 )' "
        "| perl src/popoolation2_1201/indel_filtering/identify-indel-regions.pl "
            "--input /dev/stdin "
            "--output {output.gtf} "
            "--indel-window {params.indel_window} "
            "--min-count {params.min_count}) "
        "2> {log} 1>&2"



rule sync_filter_indels:
    """Filter indels from the joint mpileup"""
    input:
        mpileups = expand(
            MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
            population=POPULATIONS,
            chromosome="{chromosome}"
        ),
        gtf = temp(SYNC_FILT + "{chromosome}.gtf")
    output:
        mpileup_fifo = temp(
            SYNC_FILT + "{chromosome}.mpileup"
        ),
        mpileup_gz = SYNC_FILT + "{chromosome}.mpileup.gz"
    params:
        mpileups_comma = "{" + ",".join(
            expand(
                MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
                population=POPULATIONS,
                chromosome="{chromosome}"
            )
        ) + "}"
    threads: 1
    log: SYNC_FILT + "{chromosome}.filter_indels.log"
    benchmark: SYNC_FILT + "{chromosome}.filter_indels.bmk"
    conda: "sync.yml"
    shell:
        "mkfifo {output.mpileup_fifo}; "
        "(cat {output.mpileup_fifo} | gzip --fast > {output.mpileup_gz} &);"
        "(eval "
            "paste <(gzip -dc {input.mpileups[0]} | cut -f 1-3) "
            "'<(gzip -dc '{params.mpileups_comma}' | cut -f 4-6 )' "
        "| perl src/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl "
            "--input /dev/stdin "
            "--gtf {input.gtf} "
            "--output {output.mpileup_fifo}) "
        "2> {log} 1>&2"


rule sync_mpileup2sync:
    """Convert joint mpileup to sync

    - mpileup2sync returns error always, and that is why there is a || true.
    - Next step requires a proper file
    """
    input:
        mpileup_gz = SYNC_FILT + "{chromosome}.mpileup.gz"
    output:
        sync = temp(SYNC_RAW + "{chromosome}.sync")  # TEMP!
    params:
        min_qual = config["popoolation2_params"]["mpileup2sync"]["min_qual"],
    threads: 1
    log: SYNC_RAW + "{chromosome}.log"
    benchmark: SYNC_RAW + "{chromosome}.json"
    resources:
        memory_gb = config["popoolation2_params"]["mpileup2sync"]["memory_gb"]
    conda: "sync.yml"
    shell:
        "(gzip --decompress --stdout {input.mpileup_gz} "
        "| java -Xmx{resources.memory_gb}g -jar src/popoolation2_1201/mpileup2sync.jar "
            "--input /dev/stdin "
            "--output {output.sync} "
            "--fastq-type sanger "
            "--min-qual {params.min_qual} "
            "--threads {threads} "
        "|| true) "
        "2> {log} 1>&2"



rule sync_subsample:
    """
    Subsample a sync file.

    Note: A proper file is required as input and output.
    """
    input:
        sync = SYNC_RAW + "{chromosome}.sync"
    output:
        sync_gz = protected(
            SYNC_SUB + "{chromosome}.sync.gz"
        ),
        sync = temp(
            SYNC_SUB + "{chromosome}.sync"
        )
    params:
        target_coverage = config["popoolation2_params"]["subsample_synchronized"]["target_coverage"],
        max_coverage = ",".join([
            config["samples"][population]["max_coverage"]
            for population in POPULATIONS
        ]),
        method =  config["popoolation2_params"]["subsample_synchronized"]["method"]
    log: SYNC_SUB + "{chromosome}.log"
    benchmark: SYNC_SUB + "{chromosome}.json"
    conda: "sync.yml"
    shell:
        "perl src/popoolation2_1201/subsample-synchronized.pl "
            "--input {input.sync} "
            "--output {output.sync} "
            "--target-coverage {params.target_coverage} "
            "--max-coverage {params.max_coverage} "
            "--method {params.method} "
        "2> {log} 1>&2 ; "
        "pigz --best --keep {output.sync} 2>> {log}"



rule sync:
    input:
        [SYNC_SUB + chromosome + ".sync.gz" for chromosome in CHROMOSOMES]
