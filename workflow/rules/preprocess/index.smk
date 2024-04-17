rule preprocess__index:
    """Index with bwa"""
    input:
        fa=REFERENCE / f"{REFERENCE_NAME}.fa.gz",
    output:
        mock=touch(PRE_INDEX / f"{REFERENCE_NAME}"),
        buckets=[
            PRE_INDEX / f"{reference_name}.{suffix}"
            for reference_name in [REFERENCE_NAME]
            for suffix in "0123 amb ann bwt.2bit.64 pac".split()
        ],
    log:
        PRE_INDEX / f"{REFERENCE_NAME}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        bwa-mem2 index \
            -p {output.mock} \
            {input.fa} \
        > {log} 2>&1
        """
