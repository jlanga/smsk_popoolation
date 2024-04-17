rule popoolation__plot__merge_values__:
    """Merge the varianve sliding results across chromosomes"""
    input:
        lambda w: [
            POP1_VS
            / f"{w.population}.{chromosome}.{w.analysis}.w{w.window}-s{w.step}.tsv.gz"
            for chromosome in CHROMOSOMES
        ],
    output:
        POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.tsv.gz",
    log:
        POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.tsv.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( zcat {input} \
        | cut -f 1,2,5 \
        | gzip \
        > {output} \
        ) 2> {log}
        """


rule popoolation__plot__merge_snps__:
    """Merge all the SNP files across chromosomes"""
    input:
        lambda w: [
            POP1_VS
            / f"{w.population}.{chromosome}.{w.analysis}.w{w.window}-s{w.step}.snps.gz"
            for chromosome in CHROMOSOMES
        ],
    output:
        POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.snps.gz",
    log:
        POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.snps.log",
    conda:
        "__environment__.yml"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


rule popoolation__plot__:
    """Plot a genome-wide result's from variance sliding"""
    input:
        tsv_gz=POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.tsv.gz",
    output:
        pdf=POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.pdf",
    log:
        POP1_PLOTS / "{population}.{analysis}.w{window}-s{step}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input} \
            --output {output.pdf} \
        2> {log}
        """


rule popoolation__plot:
    input:
        [
            POP1_PLOTS / f"{population}.{analysis}.w{window}-s{step}.pdf"
            for population in POPULATIONS
            for analysis in POP1_ANALYSES
            for window, step in POP1_WINDOW_STEP
        ],
