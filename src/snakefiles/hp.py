rule hp_table_population_chromosome:
    input:
        snps_gz =  TABLE_D + "{population}/{chromosome}.snps.gz"
    output:
        hp_gz = TABLE_HP + "{population}/{chromosome}.tsv.gz"
    params:
        
    log:
        TABLE_HP + "{population}/{chromosome}.log"
    benchmark:
        TABLE_HP + "{population}/{chromosome}.json"
    shell:
        "(pigz --decompress --stdout {input.snps_gz} "
            "| python3 src/snps_to_hp.py "
            "| pigz --best "
        "> {output.hp_gz}) "
        "2> {log}"



rule hp_plot_population:
    input:
        tsvs =expand(
            TABLE_HP + "{population}/{chromosome}.tsv.gz",
            population = "{population}",
            chromosome = CHROMOSOMES
        )
    output:
        merged_tsv_gz = PLOT_HP + "{population}.tsv.gz", 
        z_pdf = PLOT_HP + "{population}_z.pdf",
        pdf = PLOT_HP + "{population}.pdf"
    params:
        merged_tsv = PLOT_HP + "{population}.tsv"
    threads:
        1
    log:
        PLOT_HP + "{population}.log"
    benchmark:
        PLOT_HP + "{population}.json"
    shell:
        "pigz "
            "--decompress "
            "--stdout " 
            "{input.tsvs} "
        "> {params.merged_tsv} "
        "2> {log} ; "
        "Rscript src/plot_score.R "
            "none "
            "{params.merged_tsv} "
            "{output.pdf} "
        "2>> {log} ; "
        "Rscript src/plot_score.R "
            "z "
            "{params.merged_tsv} "
            "{output.z_pdf} "
        "2>> {log} ; "
        "pigz "
            "--best "
            "{params.merged_tsv} "
        "2>> {log}"