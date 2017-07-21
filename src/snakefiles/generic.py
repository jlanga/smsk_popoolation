rule index_bam:
    input: "{filename}.bam"
    output: "{filename}.bam.bai"
    shell: "samtools index {input} {output}"

rule index_cram:
    input: "{filename}.cram"
    output: "{filename}.cram.crai"
    shell: "samtools index {input} {output}"

rule index_fasta:
    input: "{filename}.fasta"
    output: "{filename}.fasta.fai"
    shell: "samtools faidx {input} {output}"

rule index_fa:
    input: "{filename}.fa"
    output: "{filename}.fa.fai"
    shell: "samtools faidx {input} {output}"
