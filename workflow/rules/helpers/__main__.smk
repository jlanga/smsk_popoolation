rule index_bam:
    """Index a bam file"""
    input:
        "{filename}.bam",
    output:
        "{filename}.bam.bai",
    log:
        "{filename}.bam.bai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input}"


rule index_cram:
    """Index a .cram file"""
    input:
        "{filename}.cram",
    output:
        "{filename}.cram.crai",
    log:
        "{filename}.cram.crai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input}"


rule index_fasta:
    """Index a .fasta"""
    input:
        "{filename}",
    output:
        "{filename}.fai",
    log:
        "{filename}.fai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools faidx {input}"
