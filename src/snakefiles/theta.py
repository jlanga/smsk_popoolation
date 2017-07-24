rule theta_table_population_chromosome:
    """
    Get the sliding Theta values.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps = temp(TABLE_THETA + "{population}/{chromosome}.snps"),
        vs = temp(TABLE_THETA + "{population}/{chromosome}.tsv"),
    params:
        mincount = config["popoolation_params"]["theta"]["mincount"],
        mincoverage = config["popoolation_params"]["theta"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["theta"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["theta"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["theta"]["poolsize"],
        stepsize = config["popoolation_params"]["theta"]["stepsize"],
        windowsize = config["popoolation_params"]["theta"]["windowsize"],
    log: TABLE_THETA + "{population}/{chromosome}.log"
    benchmark: TABLE_THETA + "{population}/{chromosome}.json"
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
            "--input <(gzip --decompress --stdout {input.mpileup_gz}) "
            "--output {output.vs} "
            "--snp-output {output.snps} "
        "2> {log} 1>&2"


rule theta_merge_vs:
    """
    Merge all tsvs into a single tsv.gz
    """
    input:
        expand(
            TABLE_THETA + "{population}/{chromosome}.tsv",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_THETA + "{population}.tsv.gz")
    threads: 8
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"


rule theta_merge_snps:
    """
    Merge all snps into a single snps.tsv
    """
    input:
        expand(
            TABLE_THETA + "{population}/{chromosome}.snps",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_THETA + "{population}.snps.gz")
    threads: 8
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"



rule theta_plot_population:
    """
    Plot a genome-wide Theta's distribution
    """
    input:
        tsv_gz = PLOT_THETA + "{population}.tsv.gz"
    output:
        z_pdf = protected(PLOT_THETA + "{population}_z.pdf"),
        pdf = protected(PLOT_THETA + "{population}.pdf")
    params:
        merged_tsv = PLOT_THETA + "{population}.tsv"
    log: PLOT_THETA + "{population}.log"
    benchmark: PLOT_THETA + "{population}.json"
    shell:
        "Rscript src/plot_score.R "
            "none "
            "{input.tsv_gz} "
            "{output.pdf} "
        "2>> {log} ; "
        "Rscript src/plot_score.R "
            "z "
            "{input.tsv_gz} "
            "{output.z_pdf} "
        "2>> {log}"
