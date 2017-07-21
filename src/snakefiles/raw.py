rule raw_make_links_pe_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward= protected(RAW + "{sample}_1.fq.gz"),
        reverse= protected(RAW + "{sample}_2.fq.gz")
    log:
        RAW + "make_links_pe_{sample}.log"
    benchmark:
        RAW + "make_links_pe_{sample}.json"
    shell:
        "ln -s $(readlink -f {input.forward}) {output.forward} 2> {log};"
        "ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}"



rule raw_make_links_se_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        single= lambda wildcards: config["samples_se"][wildcards.sample]["single"],
    output:
        single= protected(RAW + "{sample}_se.fq.gz")
    log:
        RAW + "make_links_se_{sample}.log"
    benchmark:
        RAW + "make_links_se_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.single}) {output.single} 2>  {log}
        """


rule raw_extract_genome:
    """
    Extract the fasta.gz on config.yaml into genome.fa
    """
    input:
        fa_gz = config["reference"]["dna"]
    output:
        fa = RAW + "genome.fa"
    threads:
        1
    log:
        RAW + "genome.log"
    benchmark:
        RAW + "genome.json"
    shell:
        "pigz "
            "--decompress "
            "--stdout "
            "{input.fa_gz} "
        "> {output.fa} "
        "2> {log}"
