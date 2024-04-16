rule reference__recompress__:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["references"][wildcards.genome],
    output:
        fa_gz=REFERENCE / "{genome}.fa.gz",
    log:
        REFERENCE / "{genome}.log",
    conda:
        "__environment__.yml"
    threads: 8
    shell:
        """
        ( gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference:
    input:
        REFERENCE / f"{REFERENCE_NAME}.fa.gz",
