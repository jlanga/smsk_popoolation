rule popoolation2__sync__identify_indel_regions__:
    """
    Identify indels like in mpileup_popoolation_identify_indels, but all
    together
    """
    input:
        mpileup=POP2_MPILEUP / "{chromosome}.mpileup.gz",
    output:
        gtf=temp(POP2_FILT / "{chromosome}.gtf"),
    params:
        indel_window=get_sync_indel_window,
        min_count=get_sync_indel_window,
        mpileups_comma=compose_mpileups_comma,
    log:
        POP2_FILT / "{chromosome}.identify_indels.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation2/indel_filtering/identify-indel-regions.pl \
            --input <(gzip --decompress --stdout {input.mpileup}) \
            --output {output.gtf} \
            --indel-window {params.indel_window} \
            --min-count {params.min_count} \
        2> {log} 1>&2
        """


rule popoolation2__sync__filter_indels__:
    """Filter indels from the joint mpileup
    """
    input:
        mpileup=POP2_MPILEUP / "{chromosome}.mpileup.gz",
        gtf=POP2_FILT / "{chromosome}.gtf",
    output:
        mpileup=temp(POP2_FILT / "{chromosome}.mpileup"),
        mpileup_gz=POP2_FILT / "{chromosome}.mpileup.gz",
    params:
        mpileups_comma=compose_mpileups_comma,
    log:
        POP2_FILT / "{chromosome}.filter_indels.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl  workflow/scripts/popoolation2/indel_filtering/filter-sync-by-gtf.pl \
            --input <(gzip --decompress --stdout {input.mpileup}) \
            --gtf {input.gtf} \
            --output {output.mpileup} \
        2> {log} 1>&2

        gzip --keep {output.mpileup} 2>> {log} 1>&2
        """


rule popoolation2__sync_mpileup2sync__:
    """Convert joint mpileup to sync

    - mpileup2sync returns error always, and that is why there is a || true.
    - Next step requires a proper file and uncompressed
    """
    input:
        mpileup_gz=POP2_FILT / "{chromosome}.mpileup.gz",
    output:
        sync=temp(POP2_SYNC / "{chromosome}.sync"),
    params:
        min_qual=get_sync_min_qual,
    log:
        POP2_SYNC / "{chromosome}.log",
    resources:
        memory_gb=params["popoolation2"]["subsample"]["memory_gb"],
    conda:
        "__environment__.yml"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.mpileup_gz} \
        | java -Xmx{resources.memory_gb}g -jar workflow/scripts/popoolation2/mpileup2sync.jar \
            --input /dev/stdin \
            --output {output.sync} \
            --fastq-type sanger \
            --min-qual {params.min_qual} \
            --threads {threads} \
        || true) \
        2> {log} 1>&2
        """


rule popoolation2__sync_subsample__:
    """
    Subsample a sync file.

    Note: A proper file is required as input and output.
    """
    input:
        sync=POP2_SYNC / "{chromosome}.sync",
    output:
        sync=temp(POP2_SUB / "{chromosome}.sync"),
        sync_gz=POP2_SUB / "{chromosome}.sync.gz",
    params:
        target_coverage=get_sync_target_coverage,
        max_coverage=compose_max_coverages,
        method=get_sync_subsample_method,
    log:
        POP2_SUB / "{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation2/subsample-synchronized.pl \
            --input {input.sync} \
            --output {output.sync} \
            --target-coverage {params.target_coverage} \
            --max-coverage {params.max_coverage} \
            --method {params.method} \
        2> {log} 1>&2

        gzip --keep {output.sync} 2>> {log} 1>&2
        """


rule popoolation2__sync:
    input:
        [POP2_SUB / f"{chromosome}.sync.gz" for chromosome in CHROMOSOMES],
