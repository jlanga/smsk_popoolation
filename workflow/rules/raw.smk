rule raw_extract_genome:
    """Extract the fasta.gz on config.yaml into genome.fa"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa=RAW / "genome.fa",
    log:
        RAW / "genome.log",
    conda:
        "../envs/raw.yml"
    shell:
        "gzip --decompress --stdout {input.fa_gz} > {output.fa} 2> {log}"
