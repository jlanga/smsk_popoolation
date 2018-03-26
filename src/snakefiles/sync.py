rule sync_mpileup2sync_chromosome:
    """
    Call SNPs with samtools mpileup and convert to sync.

    - mpileup2sync returns error always, and that is why there is a || true.
    - the eval thing... It could be improved someway
    """
    input:
        mpileups = expand(
            MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
            population=POPULATIONS,
            chromosome="{chromosome}"
        ),
        fa  = RAW + "genome.fa",
        fai = RAW + "genome.fa.fai",
    output:
        sync = temp(SYNC_RAW + "{chromosome}.sync")  # TEMP!
    params:
        min_qual = config["popoolation2_params"]["mpileup2sync"]["min_qual"],
        mpileups_comma = "{" + ",".join(
            expand(
                MPILEUP_RAW + "{population}/{population}.{chromosome}.mpileup.gz",
                population=POPULATIONS,
                chromosome="{chromosome}"
            )
        ) + "}",
        memory = config["popoolation2_params"]["mpileup2sync"]["memory"]
    threads: 1
    log: SYNC_RAW + "{chromosome}.log"
    benchmark: SYNC_RAW + "{chromosome}.json"
    conda: "sync.yml"
    shell:
        "(eval "
            "paste <(gzip -dc {input.mpileups[0]} | cut -f 1-3) "
            "'<(gzip -dc '{params.mpileups_comma}' | cut -f 4-6 )' "
        "| java -Xmx1G -jar src/popoolation2_1201/mpileup2sync.jar "
            "--input /dev/stdin "
            "--output {output.sync} "
            "--fastq-type sanger "
            "--min-qual {params.min_qual} "
            "--threads {threads} "
        "|| true) "
        "2> {log} 1>&2"



rule sync_subsample_chromosome:
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
        max_coverage =  config["popoolation2_params"]["subsample_synchronized"]["max_coverage"],
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
        "2> {log} ; "
        "pigz --best --keep {output.sync} 2>> {log}"



rule sync:
    input:
        [SYNC_SUB + chromosome + ".sync.gz" for chromosome in CHROMOSOMES]
