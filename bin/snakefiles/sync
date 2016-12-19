def files_for_sync_mpileup_chromosome(wildcards):
    files = [MAP_FILT + population + "/" + wildcards.chromosome + ".bam"
        for population in POPULATIONS]
    return sorted(files)

rule sync_mpileup_chromosome:
    """
    Create a mpileup file with all populations and just one chromosome
    """
    input:
        bams = files_for_sync_mpileup_chromosome,
        fa  = RAW + "genome.fa",
        fai = RAW + "genome.fa.fai"
    output:
        mpileup_gz = SYNC_MPILEUP + "{chromosome}.mpileup.gz"
    threads:
        1
    params:
        
    log:
        SYNC_MPILEUP + "{chromosome}.log"
    benchmark:
        SYNC_MPILEUP + "{chromosome}.json"
    shell:
        "(samtools mpileup "
            "-B "
            "-Q 0 "
            "-f {input.fa} "
            "{input.bams} "
            "| pigz --best ) "
        "> {output.mpileup_gz} "
        "2> {log}"



rule sync_mpileup2sync_chromosome:
    input:
        mpileup_gz = SYNC_MPILEUP + "{chromosome}.mpileup.gz"
    output:
        sync_gz = SYNC_RAW + "{chromosome}.sync.gz"
    params:
        mpileup = SYNC_MPILEUP + "{chromosome}.mpileup",
        sync = SYNC_RAW + "{chromosome}.sync",
        min_qual = config["popoolation2_params"]["mpileup2sync"]["min_qual"],    
    threads:
        8
    log:
        SYNC_RAW + "{chromosome}.log"
    benchmark:
        SYNC_RAW + "{chromosome}.json"
    shell:
        "pigz --decompress --stdout {input.mpileup_gz} > {params.mpileup} 2> {log} ; "
        "(java -jar src/popoolation2_1201/mpileup2sync.jar "
            "--input {params.mpileup} "
            "--output {params.sync} "
            "--fastq-type sanger "
            "--min-qual {params.min_qual} "
            "--threads {threads} || true ) "
        "2>> {log} ; "
        "rm {params.mpileup} ; "
        "pigz --best {params.sync} 2>> {log}"



rule sync_subsample_chromosome:
    input:
        sync_gz = SYNC_RAW + "{chromosome}.sync.gz"
    output:
        sync_gz = protected(
            SYNC_SUB + "{chromosome}.sync.gz"
        )
    params:
        sync_in = SYNC_SUB + "{chromosome}.sync.old", 
        sync_out = SYNC_SUB + "{chromosome}.sync",
        target_coverage = config["popoolation2_params"]["subsample_synchronized"]["target_coverage"],
        max_coverage =  config["popoolation2_params"]["subsample_synchronized"]["max_coverage"],
        method =  config["popoolation2_params"]["subsample_synchronized"]["method"]
    threads:
        1
    log:
        SYNC_SUB + "{chromosome}.log"
    benchmark:
        SYNC_SUB + "{chromosome}.json"
    shell:
        "pigz "
            "--decompress "
            "--stdout "
            "{input.sync_gz} "
        "> {params.sync_in} "
        "2> {log} ; "
        "perl src/popoolation2_1201/subsample-synchronized.pl "
            "--input {params.sync_in} "
            "--output {params.sync_out} "
            "--target-coverage {params.target_coverage} "
            "--max-coverage {params.max_coverage} "
            "--method {params.method} "
        "2>> {log} ; "
        "pigz --best {params.sync_out} 2>> {log} ; "
        "rm {params.sync_in} 2>> {log}"
