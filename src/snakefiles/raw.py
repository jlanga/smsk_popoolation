rule raw_make_links_pe_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.sample][wildcards.library]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.sample][wildcards.library]["reverse"]
    output:
        forward= protected(RAW + "{sample}/{library}_1.fq.gz"),
        reverse= protected(RAW + "{sample}/{library}_2.fq.gz")
    log: RAW + "{sample}/{library}/make_links_pe.log"
    benchmark: RAW + "{sample}/{library}/make_links_pe.json"
    shell:
        "ln -s $(readlink -f {input.forward}) {output.forward} 2> {log};"
        "ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}"



rule raw_make_links_se_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        single= lambda wildcards: config["samples_se"][wildcards.sample][wildcards.library]["single"],
    output:
        single= protected(RAW + "{sample}/{library}_se.fq.gz")
    log: RAW + "{sample}/{library}/make_links_se.log"
    benchmark: RAW + "{sample}_{library}/make_links_se.json"
    shell:
        "ln -s $(readlink -f {input.single}) {output.single} 2> {log}"


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
