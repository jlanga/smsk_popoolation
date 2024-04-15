rule preprocess__mpileup__:
    """Compute the mpileup and compress it"""
    input:
        cram=get_library_files_from_sample,
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        mpileup_gz=MPILEUP_RAW / "{population}" / "{population}.{chromosome}.mpileup.gz",
    log:
        MPILEUP_RAW / "{population}" /"{population}.{chromosome}.log",
    conda:
        "__environment__.yml"
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



rule preprocess__mpileup:
    input:
        [
            MPILEUP_RAW / population / f"{population}.{chromosome}.mpileup.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
        ],
