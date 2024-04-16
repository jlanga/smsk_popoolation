rule preprocess__mpileup__:
    """Compute the mpileup and compress it

    Note: don't update samtools mpileup to bcftools becaus the outputs are way different.
    """
    input:
        cram=get_library_files_from_sample,
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        mpileup_gz=PRE_MPILEUP / "{population}.{chromosome}.mpileup.gz",
    log:
        PRE_MPILEUP / "{population}.{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( samtools merge \
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
        | bgzip \
            --stdout \
        > {output.mpileup_gz} \
        ) 2> {log}
        """


rule preprocess__mpileup:
    input:
        [
            PRE_MPILEUP / f"{population}.{chromosome}.mpileup.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
        ],
