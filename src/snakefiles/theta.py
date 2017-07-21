rule theta_table_population_chromosome:
    """
    Get the Theta distribution.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps_gz = TABLE_T + "{population}/{chromosome}.snps.gz",
        vs_gz = TABLE_T + "{population}/{chromosome}.tsv.gz"
    params:
        snps = TABLE_T + "{population}/{chromosome}.snps",
        vs = TABLE_T + "{population}/{chromosome}.tsv",
        mincount = config["popoolation_params"]["theta"]["mincount"],
        mincoverage = config["popoolation_params"]["theta"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["theta"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["theta"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["theta"]["poolsize"],
        stepsize = config["popoolation_params"]["theta"]["stepsize"],
        windowsize = config["popoolation_params"]["theta"]["windowsize"],
    log: TABLE_T + "{population}/{chromosome}.log"
    benchmark: TABLE_T + "{population}/{chromosome}.json"
    shell:
        "perl src/popoolation_1.2.2/Variance-sliding.pl "
            "--measure theta "
            "--fastq-type sanger "
            "--min-count {params.mincount} "
            "--min-coverage {params.mincoverage} "
            "--max-coverage {params.maxcoverage} "
            "--min-covered-fraction {params.mincoveredfraction} "
            "--pool-size {params.poolsize} "
            "--window-size {params.windowsize} "
            "--step-size {params.stepsize} "
            "--input <(pigz --decompress --stdout {input.mpileup_gz}) "
            "--output {params.vs} "
            "--snp-output {params.snps} "
        "2> {log} 1>&2 ; "
        "pigz --best {params.snps} 2>> {log} ; "
        "pigz --best {params.vs} 2>> {log}"



rule theta_plot_population:
    input:
        tsvs =expand(
            TABLE_T + "{population}/{chromosome}.tsv.gz",
            population = "{population}",
            chromosome = CHROMOSOMES
        )
    output:
        merged_tsv_gz = PLOT_T + "{population}.tsv.gz",
        z_pdf = PLOT_T + "{population}_z.pdf",
        pdf = PLOT_T + "{population}.pdf"
    params:
        merged_tsv = PLOT_T + "{population}.tsv"
    log: PLOT_T + "{population}.log"
    benchmark: PLOT_T + "{population}.json"
    shell:
        "pigz --decompress --stdout {input.tsvs} "
            "| bash src/variance_sliding_to_genomic_score.sh "
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
        "pigz --best {params.merged_tsv} 2>> {log}"
