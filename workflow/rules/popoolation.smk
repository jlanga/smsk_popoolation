rule popoolation_variance_sliding:
    """
    Run popoolation's Variance sliding script: Tajima's D, Tajima's Theta or Pi
    """
    input:
        mpileup_gz=MPILEUP_SUB / "{population}/{population}.{chromosome}.mpileup.gz",
    output:
        snps=POPOOLATION_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.snps",
        vs=POPOOLATION_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.tsv",
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
        POPOOLATION_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.log",
    conda:
        "../envs/popoolation.yml"
    shell:
        """
        perl workflow/scripts/popoolation_1.2.2/Variance-sliding.pl \
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


rule popoolation_merge_variance_sliding:
    """Merge results across chromosomes"""
    input:
        expand(
            POPOOLATION_TABLES / "{analysis}/{population}.{chromosome}.{analysis}.tsv",
            chromosome=CHROMOSOMES,
            population=["{population}"],
            analysis=["{analysis}"],
        ),
    output:
        protected(POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.tsv.gz"),
    log:
        POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.merge_vs.log",
    threads: 24
    conda:
        "../envs/popoolation.yml"
    shell:
        "(bash workflow/scripts/variance_sliding_to_genomic_score.sh {input} "
        "| pigz --best --processes {threads} > {output}) "
        "2> {log}"


rule popoolation_merge_snps:
    """Merge snps"""
    input:
        expand(
            POPOOLATION_TABLES
            / "{analysis}/{population}.{chromosome}.{analysis}.snps",
            chromosome=CHROMOSOMES,
            population=["{population}"],
            analysis=["{analysis}"],
        ),
    output:
        protected(POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.snps.gz"),
    threads: 8
    log:
        POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.merge_snps.log",
    conda:
        "../envs/popoolation.yml"
    shell:
        "pigz --processes {threads} --best --stdout {input} > {output}"


rule popoolation_plots:
    """Plot a genome-wide result's from variance sliding"""
    input:
        tsv_gz=POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.tsv.gz",
    output:
        pdf=protected(POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.pdf"),
    log:
        POPOOLATION_PLOTS / "{analysis}/{population}.{analysis}.plot.log",
    conda:
        "../envs/popoolation.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input.tsv_gz} \
            --output {output.pdf} \
        2> {log}
        """


rule popoolation_d:
    input:
        expand(POPOOLATION_PLOTS / "D/{population}.D.pdf", population=POPULATIONS),


rule popoolation_pi:
    input:
        expand(POPOOLATION_PLOTS / "pi/{population}.pi.pdf", population=POPULATIONS),


rule popoolation_theta:
    input:
        expand(
            POPOOLATION_PLOTS / "theta/{population}.theta.pdf", population=POPULATIONS
        ),


rule popoolation:
    input:
        rules.popoolation_d.input,
        rules.popoolation_pi.input,
        rules.popoolation_theta.input,
