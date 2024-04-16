rule popoolation__variance_sliding__:
    """
    Run popoolation's Variance sliding script: Tajima's D, Tajima's Theta or Pi
    """
    input:
        mpileup_gz=POP1_SUB / "{population}/{population}.{chromosome}.mpileup.gz",
    output:
        snps=temp(POP1_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.snps"),
        vs=temp(POP1_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.tsv"),
    params:
        measure="{analysis}",
        min_count=get_popoolation_min_count,
        min_coverage=get_popoolation_min_coverage,
        max_coverage=get_popoolation_max_coverage,
        min_covered_fraction=get_popoolation_min_covered_fraction,
        pool_size=get_pool_size,
        step_size=get_popoolation_step_size,
        window_size=get_popoolation_window_size,
    log:
        POP1_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation/Variance-sliding.pl \
            --measure {params.measure} \
            --fastq-type sanger \
            --min-count {params.min_count} \
            --min-coverage {params.min_coverage} \
            --max-coverage {params.max_coverage} \
            --min-covered-fraction {params.min_covered_fraction} \
            --pool-size {params.pool_size} \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --input <(gzip --decompress --stdout {input.mpileup_gz}) \
            --output {output.vs} \
            --snp-output {output.snps} \
        2> {log} 1>&2
        """


rule popoolation__variance_sliding__merge_values__:
    """Merge results across chromosomes"""
    input:
        expand(
            POP1_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.tsv",
            chromosome=CHROMOSOMES,
            population=["{population}"],
            analysis=["{analysis}"],
        ),
    output:
        POP1_PLOTS / "{analysis}/{population}.{analysis}.tsv.gz",
    log:
        POP1_PLOTS / "{analysis}/{population}.{analysis}.merge_vs.log",
    threads: 24
    conda:
        "__environment__.yml"
    shell:
        """
        ( bash workflow/scripts/variance_sliding_to_genomic_score.sh {input} \
        | pigz \
            --best \
            --processes {threads} \
        > {output} \
        ) 2> {log}
        """


rule popoolation__variance_sliding__merge_values:
    input:
        [
            POP1_PLOTS / analysis / f"{population}.{analysis}.tsv.gz"
            for analysis in ["D", "pi", "theta"]
            for population in POPULATIONS
        ],


rule popoolation__variance_sliding__merge_snps__:
    """Merge snps"""
    input:
        expand(
            POP1_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.snps",
            chromosome=CHROMOSOMES,
            population=["{population}"],
            analysis=["{analysis}"],
        ),
    output:
        POP1_PLOTS / "{analysis}/{population}.{analysis}.snps.gz",
    threads: 8
    log:
        POP1_PLOTS / "{analysis}/{population}.{analysis}.merge_snps.log",
    conda:
        "__environment__.yml"
    shell:
        "pigz --processes {threads} --best --stdout {input} > {output}"


rule popoolation__variance_sliding__merge_snps:
    input:
        [
            POP1_PLOTS / analysis / f"{population}.{analysis}.snps.gz"
            for analysis in ["D", "pi", "theta"]
            for population in POPULATIONS
        ],


rule popoolation__variance_sliding:
    input:
        rules.popoolation__variance_sliding__merge_values.input,
        rules.popoolation__variance_sliding__merge_snps.input,
