rule popoolation__plot__:
    """Plot a genome-wide result's from variance sliding"""
    input:
        tsv_gz=POP1_PLOTS / "{analysis}/{population}.{analysis}.tsv.gz",
    output:
        pdf=POP1_PLOTS / "{analysis}/{population}.{analysis}.pdf",
    log:
        POP1_PLOTS / "{analysis}/{population}.{analysis}.plot.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input.tsv_gz} \
            --output {output.pdf} \
        2> {log}
        """


rule popoolation__plot:
    input:
        expand(
            POP1_PLOTS / "{analysis}/{population}.{analysis}.pdf",
            analysis=["D", "pi", "theta"],
            population=POPULATIONS,
        ),