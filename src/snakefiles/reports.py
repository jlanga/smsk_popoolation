rule reports_fastqc:
    input: "{filename}.fq.gz"
    output: "{filename}_fastqc.html", "{filename}_fastqc.zip"
    log: "{filename}.fastqc.log"
    benchmark: "{filename}.fastqc.bmk"
    conda: "reports.yml"
    shell: "fastqc --noextract --nogroup {input} 2> {log} 1>&2"


rule reports_samtools_stats:
    input: "{filename}.cram"
    output: "{filename}.stats.tsv"
    log: "{filename}.stats.log"
    benchmark: "{filename}.stats.bmk"
    conda: "reports.yml"
    shell: "samtools stats {input} > {output} 2>&1"


rule reports_samtools_flagstat:
    input: "{filename}.cram"
    output: "{filename}.flagstat.txt"
    log: "{filename}.flagstat.log"
    benchmark: "{filename}.flagstat.bmk"
    conda: "reports.yml"
    shell: "samtools flagstat {input} > {output} 2> {log}"


rule reports_samtools_idxstats:
    input:
        cram = "{filename}.cram",
        crai = "{filename}.cram.crai"
    output: "{filename}.idxstats.txt"
    log: "{filename}.idxstats.log"
    benchmark: "{filename}.idxstats.bmk"
    conda: "reports.yml"
    shell: "samtools idxstats {input} > {output} 2> {log}"


rule reports:
    input:
        expand(
            MAP_RAW + population + "." + library + "." + analysis
            for population, library in (
                samples
                [["population", "library"]]
                .values.tolist()
            )
            for analysis in "stats.tsv flagstat.txt idxstats.txt".split(" ")
        )
