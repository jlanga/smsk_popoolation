rule reports_fastqc:
    input:
        "{filename}.fq.gz",
    output:
        "{filename}_fastqc.html",
        "{filename}_fastqc.zip",
    log:
        "{filename}.fastqc.log",
    conda:
        "../envs/reports.yml"
    shell:
        "fastqc --noextract --nogroup {input} 2> {log} 1>&2"


rule reports_samtools_stats:
    input:
        "{filename}.cram",
    output:
        "{filename}.stats.tsv",
    log:
        "{filename}.stats.log",
    conda:
        "../envs/reports.yml"
    shell:
        "samtools stats {input} > {output} 2>&1"


rule reports_samtools_flagstat:
    input:
        "{filename}.cram",
    output:
        "{filename}.flagstat.txt",
    log:
        "{filename}.flagstat.log",
    conda:
        "../envs/reports.yml"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"


rule reports_samtools_idxstats:
    input:
        cram="{filename}.cram",
        crai="{filename}.cram.crai",
    output:
        "{filename}.idxstats.txt",
    log:
        "{filename}.idxstats.log",
    conda:
        "../envs/reports.yml"
    shell:
        "samtools idxstats {input.cram} > {output} 2> {log}"


rule reports:
    input:
        expand(
            MAP_RAW / f"{population}.{library}.{analysis}"
            for population, library in POPULATION_LIBRARY
            for analysis in "stats.tsv flagstat.txt idxstats.txt".split(" ")
        ),
