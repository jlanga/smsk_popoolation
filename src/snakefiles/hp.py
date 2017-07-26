rule hp_table_population_chromosome:
    """
    Build the hp table for a population in a chromosome.
    """
    input:
        snps_gz =  PLOT_D + "{population}.snps.gz"
    output:
        hp_gz = protected(TABLE_HP + "{population}.tsv.gz")
    log: TABLE_HP + "{population}.log"
    benchmark: TABLE_HP + "{population}.json"
    shell:
        "(gzip --decompress --stdout {input.snps_gz} "
        "| python3 src/snps_to_hp.py "
        "| pigz --best "
        "> {output.hp_gz}) "
        "2> {log}"



rule hp_plot_population:
    """
    Plot the genome-wide H_p distribution of a population
    """
    input:
        hp_gz = TABLE_HP + "{population}.tsv.gz"
    output:
        z_pdf = PLOT_HP + "{population}_z.pdf",
        pdf = PLOT_HP + "{population}.pdf"
    threads:
        1
    log: PLOT_HP + "{population}.log"
    benchmark: PLOT_HP + "{population}.json"
    shell:
        "Rscript src/plot_score.R "
            "--input {input.hp_gz} "
            "--output {output.pdf} "
        "2>> {log} ; "
        "Rscript src/plot_score.R "
            "--normalize "
            "--input {input.hp_gz} "
            "--output {output.z_pdf} "
        "2>> {log}"
