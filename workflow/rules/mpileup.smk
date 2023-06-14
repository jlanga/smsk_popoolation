rule mpileup_convert:
    """Compute the mpileup and compress it"""
    input:
        cram=get_library_files_from_sample,
        fa=RAW / "genome.fa",
        fai=RAW / "genome.fa.fai",
    output:
        mpileup_gz=MPILEUP_RAW / "{population}/{population}.{chromosome}.mpileup.gz",
    log:
        MPILEUP_RAW / "{population}/{population}.{chromosome}.log",
    conda:
        "../envs/mpileup.yml"
    shell:
        """
        (samtools merge \
            -u \
            --reference {input.fa} \
            - \
            {input.cram} \
        | samtools mpileup \
            -a \
            --no-BAQ \
            --min-BQ 0 \
            --fasta-ref {input.fa} \
            - \
        | gzip \
            --best \
        > {output.mpileup_gz} \
        ) 2> {log}
        """


rule mpileup_popoolation_identify_indels:
    """
    Get a GTF with the indels present.
    """
    input:
        mpileup_gz=MPILEUP_RAW / "{population}/{population}.{chromosome}.mpileup.gz",
    output:
        gtf=temp(MPILEUP_FILT / "{population}/{population}.{chromosome}.gtf"),
    params:
        indel_window=get_indel_window,
        min_count=get_indel_min_count,
    log:
        MPILEUP_FILT / "{population}/{population}.{chromosome}.gtf.log",
    conda:
        "../envs/mpileup.yml"
    shell:
        """
        perl workflow/scripts/popoolation_1.2.2/basic-pipeline/\
identify-genomic-indel-regions.pl \
            --input <(gzip --decompress --stdout {input.mpileup_gz}) \
            --output {output.gtf} \
            --indel-window {params.indel_window} \
            --min-count {params.min_count} \
        2>{log} 1>&2
        """


rule mpileup_popoolation_filter_indels:
    """
    Filter indels from an mpileup given a GTF files with their coordinates.
    Compress results.

    Both fifo and fifo.gz will be deleted.
    """
    input:
        mpileup_gz=MPILEUP_RAW / "{population}/{population}.{chromosome}.mpileup.gz",
        gtf=MPILEUP_FILT / "{population}/{population}.{chromosome}.gtf",
    output:
        mpileup_fifo=temp(
            MPILEUP_FILT / "{population}/{population}.{chromosome}.mpileup"
        ),
        mpileup_gz=temp(
            MPILEUP_FILT / "{population}/{population}.{chromosome}.mpileup.gz"
        ),
    log:
        MPILEUP_FILT / "{population}/{population}.{chromosome}.mpileup.log",
    conda:
        "../envs/mpileup.yml"
    shell:
        """
        mkfifo {output.mpileup_fifo}

        (cat {output.mpileup_fifo} | gzip --fast > {output.mpileup_gz} &)

        perl workflow/scripts/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
            --input <(gzip --decompress --stdout {input.mpileup_gz}) \
            --gtf {input.gtf} \
            --output {output.mpileup_fifo} \
        2> {log} 1>&2
        """


rule mpileup_popoolation_subsample:
    """
    Perform the subsampling step. Compress results as they are generated
    through a FIFO
    """
    input:
        mpileup=MPILEUP_FILT / "{population}/{population}.{chromosome}.mpileup.gz",
    output:
        mpileup_fifo=temp(
            MPILEUP_SUB / "{population}/{population}.{chromosome}.mpileup"
        ),
        mpileup_gz=protected(
            MPILEUP_SUB / "{population}/{population}.{chromosome}.mpileup.gz"
        ),
    params:
        min_qual=get_subsample_min_qual,
        method=get_subsample_method,
        max_coverage=get_subsample_max_coverage,
        target_coverage=get_subsample_target_coverage,
    log:
        MPILEUP_SUB / "{population}/{population}.{chromosome}.log",
    conda:
        "../envs/mpileup.yml"
    shell:
        """
        rm -f {output.mpileup_fifo}
        mkfifo {output.mpileup_fifo}

        (cat {output.mpileup_fifo} | gzip --best > {output.mpileup_gz} &)

        perl workflow/scripts/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl \
            --min-qual {params.min_qual} \
            --method {params.method} \
            --max-coverage {params.max_coverage} \
            --fastq-type sanger \
            --target-coverage {params.target_coverage} \
            --input <(gzip -dc {input.mpileup}) \
            --output {output.mpileup_fifo} \
        2> {log} 1>&2
        """


rule mpileup:
    input:
        [
            MPILEUP_SUB / f"{population}/{population}.{chromosome}.mpileup.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
        ],
