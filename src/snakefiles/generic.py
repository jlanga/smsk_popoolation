rule index_bam:
    """Index a bam file"""
    input: "{filename}.bam"
    output: "{filename}.bam.bai"
    shell: "samtools index {input}"

rule index_cram:
    """Index a .cram file"""
    input: "{filename}.cram"
    output: "{filename}.cram.crai"
    shell: "samtools index {input}"

rule index_fasta:
    """Index a .fasta"""
    input: "{filename}.fasta"
    output: "{filename}.fasta.fai"
    shell: "samtools faidx {input}"

rule index_fa:
    """Index a .fa"""
    input: "{filename}.fa"
    output: "{filename}.fa.fai"
    shell: "samtools faidx {input}"

rule decompress_mpileup:
    """Decompress temporarily a .mpileup.gz"""
    input: "{filename}.mpileup.gz"
    output: temp("{filename}.mpileup")
    shell: "pigz --decompress --keep {input}"

rule decompress_tsv:
    """Decompress temporarily a .tsv"""
    input: "{filename}.tsv.gz"
    output: temp("{filename}.tsv")
    shell: "pigz --decompress --keep {input}"

rule decompress_sync:
    """Decompress temporarily a .sync"""
    input: "{filename}.sync.gz"
    output: "{filename}.sync"
    shell: "pigz --decompress --keep {input}"

# rule join_compressed_tsvs:
#     """Join multiple tsvs into one, and compress it"""
#     input: expand(
#         "{prefix}/{chromosome}.tsv.gz",
#         chromosome = CHROMOSOMES
#     )
#     output: "{prefix}.tsv.gz"
#     shell:
#         "pigz --decompress --stdout {input} "
#         "| pigz --best > {output}"
