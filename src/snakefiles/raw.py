rule raw_make_links_pe_sample:
    """
    Make a link to the original file, with a prettier name than default.
    """
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.sample][wildcards.library]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.sample][wildcards.library]["reverse"]
    output:

        forward= RAW + "{sample}/{library}_1.fq.gz",
        reverse= RAW + "{sample}/{library}_2.fq.gz"
    shell:
        "ln --symbolic $(readlink --canonicalize {input.forward}) {output.forward}; "
        "ln --symbolic $(readlink --canonicalize {input.reverse}) {output.reverse}"



rule raw_make_links_se_sample:
    """
    Make a link to the original file, with a prettier name than default.
    """
    input:
        single= lambda wildcards: config["samples_se"][wildcards.sample][wildcards.library]["single"],
    output:
        single= RAW + "{sample}/{library}_se.fq.gz"
    shell:
        "ln --symbolic $(readlink --canonicalize {input.single}) {output.single}"


rule raw_extract_genome:
    """
    Extract the fasta.gz on config.yaml into genome.fa
    """
    input:
        fa_gz = config["reference"]["dna"]
    output:
        fa = RAW + "genome.fa"
    log: RAW + "genome.log"
    benchmark: RAW + "genome.json"
    shell:
        "pigz "
            "--decompress "
            "--stdout "
            "{input.fa_gz} "
        "> {output.fa} "
        "2> {log}"
