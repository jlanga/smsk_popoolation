rule preprocess__mpileup__:
    """Compute the mpileup and compress it

    Join multiple libraries
    Don't update samtools to bcftools:
    - samtools mpileup produces mpileup format
    - bcftools mpileup produces VCF format
    """
    input:
        bams=get_library_files_from_sample,
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
        fai=REFERENCE / f"{REFERENCE_NAME}.fa.gz.fai",
    output:
        mpileup_gz=PRE_MPILEUP / "{population}" / "{chromosome}.mpileup.gz",
    log:
        PRE_MPILEUP / "{population}" / "{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( samtools merge \
            -u \
            --reference {input.fa} \
            - \
            {input.bams} \
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
            PRE_MPILEUP / population / f"{chromosome}.mpileup.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
        ],
