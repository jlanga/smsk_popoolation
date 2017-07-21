rule tajimad_table_population_chromosome:
    """
    Get the sliding Tajima's D values.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps_gz = TABLE_D + "{population}/{chromosome}.snps.gz",
        vs_gz = TABLE_D + "{population}/{chromosome}.tsv.gz"
    params:
        snps = TABLE_D + "{population}/{chromosome}.snps",
        vs = TABLE_D + "{population}/{chromosome}.tsv",
        mincount = config["popoolation_params"]["tajimad"]["mincount"],
        mincoverage = config["popoolation_params"]["tajimad"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["tajimad"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["tajimad"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["tajimad"]["poolsize"],
        stepsize = config["popoolation_params"]["tajimad"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimad"]["windowsize"],
    log: TABLE_D + "{population}/{chromosome}.log"
    benchmark: TABLE_D + "{population}/{chromosome}.json"
    shell:
        "perl src/popoolation_1.2.2/Variance-sliding.pl "
            "--measure D "
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



rule tajimad_plot_population:
    """
    Plot a genome-wide Tajima's D distribution
    """
    input:
        tsvs =expand(
            TABLE_D + "{population}/{chromosome}.tsv.gz",
            population = "{population}",
            chromosome = CHROMOSOMES
        )
    output:
        merged_tsv_gz = PLOT_D + "{population}.tsv.gz",
        z_pdf = PLOT_D + "{population}_z.pdf",
        pdf = PLOT_D + "{population}.pdf"
    params:
        merged_tsv = PLOT_D + "{population}.tsv"
    log: PLOT_D + "{population}.log"
    benchmark: PLOT_D + "{population}.json"
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
