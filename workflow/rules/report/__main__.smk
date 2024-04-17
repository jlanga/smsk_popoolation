rule report__fastqc__:
    input:
        "{filename}.fq.gz",
    output:
        "{filename}_fastqc.html",
        "{filename}_fastqc.zip",
    log:
        "{filename}.fastqc.log",
    conda:
        "__environment__.yml"
    shell:
        "fastqc --noextract --nogroup {input} 2> {log} 1>&2"


rule report__samtools_stats__:
    input:
        "{filename}.cram",
    output:
        "{filename}.stats.tsv",
    log:
        "{filename}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats {input} > {output} 2>&1"


rule report__samtools_flagstat__:
    input:
        "{filename}.cram",
    output:
        "{filename}.flagstat.txt",
    log:
        "{filename}.flagstat.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"


rule report__samtools_idxstats__:
    input:
        cram="{filename}.cram",
        crai="{filename}.cram.crai",
    output:
        "{filename}.idxstats.txt",
    log:
        "{filename}.idxstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools idxstats {input.cram} > {output} 2> {log}"


rule report:
    input:
        [
            PRE_MAP / f"{population}.{library}.{analysis}"
            for population, library in POPULATION_LIBRARY
            for analysis in "stats.tsv flagstat.txt idxstats.txt".split(" ")
        ],
