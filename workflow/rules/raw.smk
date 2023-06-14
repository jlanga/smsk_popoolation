rule raw_make_links_pe:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=RAW / "{population}.{library}_1.fq.gz",
        reverse_=RAW / "{population}.{library}_2.fq.gz",
    log:
        RAW / "{population}.{library}.log",
    conda:
        "../envs/raw.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


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
