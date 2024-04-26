rule preprocess__index:
    """Index with bwa"""
    input:
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        buckets=multiext(
            str(PRE_INDEX / f"{REFERENCE_NAME}."),
            "0123",
            "amb",
            "ann",
            "bwt.2bit.64",
            "pac",
        ),
    log:
        PRE_INDEX / f"{REFERENCE_NAME}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=PRE_INDEX / f"{REFERENCE_NAME}",
    shell:
        """
        bwa-mem2 index \
            -p {params.prefix} \
            {input.fa} \
        > {log} 2>&1
        """
