def get_reads(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    forward_, reverse_ = samples[
        (samples["population"] == pop) & (samples["library"] == lib)
    ][["forward", "reverse"]].values.tolist()[0]
    return forward_, reverse_


rule raw_make_links_pe:
    """Make a link to the original file, with a prettier name than default"""
    input:
        get_reads,
    output:
        forward_=RAW + "{population}.{library}_1.fq.gz",
        reverse_=RAW + "{population}.{library}_2.fq.gz",
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input[0]}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input[1]}) {output.reverse_}
        """


rule raw_extract_genome:
    """Extract the fasta.gz on config.yaml into genome.fa"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa=RAW + "genome.fa",
    log:
        RAW + "genome.log",
    benchmark:
        RAW + "genome.json"
    shell:
        "gzip --decompress --stdout {input.fa_gz} > {output.fa} 2> {log}"
