rule popoolation__hp__compute__:
    """
    Build the hp table for a population in a chromosome.
    """
    input:
        snps_gz=POP1_PLOTS / "{population}.D.w{window}-s{step}.snps.gz",
    output:
        hp_gz=POP1_HP / "{population}.w{window}-s{step}.tsv.gz",
    log:
        POP1_HP / "{population}.w{window}-s{step}.log",
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
        hp_gz=POP1_HP / "{population}.w{window}-s{step}.tsv.gz",
    output:
        pdf=POP1_HP / "{population}.w{window}-s{step}.pdf",
    log:
        POP1_HP / "{population}.w{window}-s{step}.log",
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
        [
            POP1_HP / f"{population}.w{window}-s{step}.pdf"
            for population in POPULATIONS
            for window, step in POP1_WINDOW_STEP
        ],
