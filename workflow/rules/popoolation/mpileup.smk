rule popoolation__mpileup__identify_indels__:
    """
    Get a GTF with the indels present.
    """
    input:
        mpileup_gz=PRE_MPILEUP / "{population}.{chromosome}.mpileup.gz",
    output:
        gtf=temp(POP1_FILT / "{population}.{chromosome}.gtf"),
    params:
        indel_window=get_indel_window,
        min_count=get_indel_min_count,
    log:
        POP1_FILT / "{population}.{chromosome}.gtf.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation/basic-pipeline/\
identify-genomic-indel-regions.pl \
            --input <(pigz --decompress --stdout {input.mpileup_gz}) \
            --output {output.gtf} \
            --indel-window {params.indel_window} \
            --min-count {params.min_count} \
        2>{log} 1>&2
        """


rule popoolation__mpileup__filter_indels__:
    """
    Filter indels from an mpileup given a GTF files with their coordinates.
    Compress results.

    Both fifo and fifo.gz will be deleted.
    """
    input:
        mpileup_gz=PRE_MPILEUP / "{population}.{chromosome}.mpileup.gz",
        gtf=POP1_FILT / "{population}.{chromosome}.gtf",
    output:
        mpileup_fifo=temp(POP1_FILT / "{population}.{chromosome}.mpileup"),
        mpileup_gz=temp(POP1_FILT / "{population}.{chromosome}.mpileup.gz"),
    log:
        POP1_FILT / "{population}.{chromosome}.mpileup.log",
    conda:
        "__environment__.yml"
    shell:
        """
        mkfifo {output.mpileup_fifo}

        (cat {output.mpileup_fifo} | pigz --fast > {output.mpileup_gz} &)

        perl workflow/scripts/popoolation/basic-pipeline/filter-pileup-by-gtf.pl \
            --input <(pigz --decompress --stdout {input.mpileup_gz}) \
            --gtf {input.gtf} \
            --output {output.mpileup_fifo} \
        2> {log} 1>&2
        """


rule popoolation__mpileup__subsample__:
    """
    Perform the subsampling step. Compress results as they are generated
    through a FIFO
    """
    input:
        mpileup=POP1_FILT / "{population}.{chromosome}.mpileup.gz",
    output:
        mpileup_fifo=temp(POP1_SUB / "{population}.{chromosome}.mpileup"),
        mpileup_gz=POP1_SUB / "{population}.{chromosome}.mpileup.gz",
    params:
        min_qual=get_subsample_min_qual,
        method=get_subsample_method,
        max_coverage=get_subsample_max_coverage,
        target_coverage=get_subsample_target_coverage,
    log:
        POP1_SUB / "{population}.{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        rm -f {output.mpileup_fifo}
        mkfifo {output.mpileup_fifo}

        pigz --fast --stdout {output.mpileup_fifo} > {output.mpileup_gz} & 2> {log}

        perl workflow/scripts/popoolation/basic-pipeline/subsample-pileup.pl \
            --min-qual {params.min_qual} \
            --method {params.method} \
            --max-coverage {params.max_coverage} \
            --fastq-type sanger \
            --target-coverage {params.target_coverage} \
            --input <(pigz --decompress --stdout {input.mpileup}) \
            --output {output.mpileup_fifo} \
        2>> {log} 1>&2
        """


rule popoolation__mpileup:
    input:
        [
            POP1_SUB / f"{population}.{chromosome}.mpileup.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
        ],
