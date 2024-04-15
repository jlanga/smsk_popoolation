rule sync_identify_indels:
    """
    Identify indels like in mpileup_popoolation_identify_indels, but all
    together
    """
    input:
        mpileups=expand(
            MPILEUP_RAW / "{population}/{population}.{chromosome}.mpileup.gz",
            population=POPULATIONS,
            chromosome="{chromosome}",
        ),
    output:
        gtf=temp(SYNC_FILT / "{chromosome}.gtf"),
    params:
        indel_window=get_sync_indel_window,
        min_count=get_sync_indel_window,
        mpileups_comma=compose_mpileups_comma,
    log:
        SYNC_FILT / "{chromosome}.identify_indels.log",
    conda:
        "../envs/sync.yml"
    shell:
        """
        (eval \
            paste <(gzip -dc {input.mpileups[0]} | cut -f 1-3) \
            '<(gzip -dc '{params.mpileups_comma}' | cut -f 4-6 )' \
        | perl workflow/scripts/popoolation2_1201/indel_filtering/\
identify-indel-regions.pl \
            --input /dev/stdin \
            --output {output.gtf} \
            --indel-window {params.indel_window} \
            --min-count {params.min_count}) \
        2> {log} 1>&2
        """


rule sync_filter_indels:
    """Filter indels from the joint mpileup"""
    input:
        mpileups=expand(
            MPILEUP_RAW / "{population}/{population}.{chromosome}.mpileup.gz",
            population=POPULATIONS,
            chromosome="{chromosome}",
        ),
        gtf=SYNC_FILT / "{chromosome}.gtf",
    output:
        mpileup_fifo=temp(SYNC_FILT / "{chromosome}.mpileup"),
        mpileup_gz=SYNC_FILT / "{chromosome}.mpileup.gz",
    params:
        mpileups_comma=compose_mpileups_comma,
    log:
        SYNC_FILT / "{chromosome}.filter_indels.log",
    conda:
        "../envs/sync.yml"
    shell:
        """
        mkfifo {output.mpileup_fifo}

        (cat {output.mpileup_fifo} | gzip --fast > {output.mpileup_gz} &)

        (eval \
            paste <(gzip -dc {input.mpileups[0]} | cut -f 1-3) \
            '<(gzip -dc '{params.mpileups_comma}' | cut -f 4-6 )' \
        | perl  workflow/scripts/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl \
            --input /dev/stdin \
            --gtf {input.gtf} \
            --output {output.mpileup_fifo}) \
        2> {log} 1>&2
        """


rule sync_mpileup2sync:
    """Convert joint mpileup to sync

    - mpileup2sync returns error always, and that is why there is a || true.
    - Next step requires a proper file
    """
    input:
        mpileup_gz=SYNC_FILT / "{chromosome}.mpileup.gz",
    output:
        sync=temp(SYNC_MPILEUP2SYNC / "{chromosome}.sync"),  # TEMP!
    params:
        min_qual=get_sync_min_qual,
    # Required for java not doing more than needed
    threads: 1
    log:
        SYNC_MPILEUP2SYNC / "{chromosome}.log",
    resources:
        memory_gb=params["popoolation2"]["subsample"]["memory_gb"],
    conda:
        "../envs/sync.yml"
    shell:
        """
        (gzip --decompress --stdout {input.mpileup_gz} \
        | java -Xmx{resources.memory_gb}g -jar workflow/scripts/popoolation2_1201/mpileup2sync.jar \
            --input /dev/stdin \
            --output {output.sync} \
            --fastq-type sanger \
            --min-qual {params.min_qual} \
            --threads {threads} \
        || true) \
        2> {log} 1>&2
        """


rule sync_subsample:
    """
    Subsample a sync file.

    Note: A proper file is required as input and output.
    """
    input:
        sync=SYNC_MPILEUP2SYNC / "{chromosome}.sync",
    output:
        sync=temp(SYNC_SUBSAMPLED / "{chromosome}.sync"),
    params:
        target_coverage=get_sync_target_coverage,
        max_coverage=compose_max_coverages,
        method=get_sync_subsample_method,
    log:
        SYNC_SUBSAMPLED / "{chromosome}.log",
    conda:
        "../envs/sync.yml"
    shell:
        """
        perl workflow/scripts/popoolation2_1201/subsample-synchronized.pl \
            --input {input.sync} \
            --output {output.sync} \
            --target-coverage {params.target_coverage} \
            --max-coverage {params.max_coverage} \
            --method {params.method} \
        2> {log} 1>&2
        """


rule sync_compress:
    input:
        sync=SYNC_SUBSAMPLED / "{chromosome}.sync",
    output:
        sync_gz=SYNC_SUBSAMPLED / "{chromosome}.sync.gz",
    log:
        SYNC_SUBSAMPLED / "{chromosome}.compressed.log",
    conda:
        "../envs/sync.yml"
    threads: 4
    shell:
        "pigz --best --keep {output.sync}"


rule sync:
    input:
        [SYNC_SUBSAMPLED / f"{chromosome}.sync.gz" for chromosome in CHROMOSOMES],
