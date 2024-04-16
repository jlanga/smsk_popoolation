rule popoolation__hp__compute__:
    """
    Build the hp table for a population in a chromosome.
    """
    input:
        snps_gz=POP1_TABLES / "{population}.D.snps.gz",
    output:
        hp_gz=HP_TABLES / "{population}.tsv.gz",
    log:
        HP_TABLES / "{population}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.snps_gz} \
        | python3 workflow/scripts/snps_to_hp.py \
        | pigz \
            --stdout \
        > {output.hp_gz} \
        ) 2> {log}
        """


rule popoolation__hp__plot__:
    """
    Plot the genome-wide H_p distribution of a population
    """
    input:
        hp_gz=HP_TABLES / "{population}.tsv.gz",
    output:
        pdf=HP_PLOTS / "{population}.pdf",
    log:
        HP_PLOTS / "{population}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input.hp_gz} \
            --output {output.pdf} \
            --logarithm \
        2> {log}
        """


rule popoolation__hp:
    input:
        [HP_PLOTS / f"{population}.pdf" for population in POPULATIONS],
