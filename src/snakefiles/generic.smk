rule index_bam:
    """Index a bam file"""
    input: "{filename}.bam"
    output: "{filename}.bam.bai"
    conda: "generic.yml"
    shell: "samtools index {input}"

rule index_cram:
    """Index a .cram file"""
    input: "{filename}.cram"
    output: "{filename}.cram.crai"
    conda: "generic.yml"
    shell: "samtools index {input}"

rule index_fasta:
    """Index a .fasta"""
    input: "{filename}"
    output: "{filename}.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"
