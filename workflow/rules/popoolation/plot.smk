rule popoolation__plot__merge_values__:
    """Merge the varianve sliding results across chromosomes"""
    input:
        lambda w: [
            POP1_VS / f"{population}.{chromosome}.{analysis}.tsv.gz"
            for population in [w.population]
            for chromosome in CHROMOSOMES
            for analysis in [w.analysis]
        ],
    output:
        POP1_PLOTS / "{population}.{analysis}.tsv.gz",
    log:
        POP1_PLOTS / "{population}.{analysis}.merge_vs.log",
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
            POP1_VS / f"{population}.{chromosome}.{analysis}.snps.gz"
            for population in [w.population]
            for chromosome in CHROMOSOMES
            for analysis in [w.analysis]
        ],
    output:
        POP1_PLOTS / "{population}.{analysis}.snps.gz",
    log:
        POP1_PLOTS / "{population}.{analysis}.merge_snps.log",
    conda:
        "__environment__.yml"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


rule popoolation__plot__:
    """Plot a genome-wide result's from variance sliding"""
    input:
        tsv_gz=POP1_PLOTS / "{population}.{analysis}.tsv.gz",
    output:
        pdf=POP1_PLOTS / "{population}.{analysis}.pdf",
    log:
        POP1_PLOTS / "{population}.{analysis}.plot.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input} \
            --output {output.pdf} \
        2> {log}
        """


rule popoolation__plot__d:
    input:
        [POP1_PLOTS / f"{population}.D.pdf" for population in POPULATIONS],


rule popoolation__plot__pi:
    input:
        [POP1_PLOTS / f"{population}.pi.pdf" for population in POPULATIONS],


rule popoolation__plot__theta:
    input:
        [POP1_PLOTS / f"{population}.theta.pdf" for population in POPULATIONS],


rule popoolation__plot:
    input:
        rules.popoolation__plot__d.input,
        rules.popoolation__plot__pi.input,
        rules.popoolation__plot__theta.input,
