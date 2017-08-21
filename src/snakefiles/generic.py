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
