rule tajimapi_table_population_chromosome:
    """
    Get the sliding Tajima's D values.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps = temp(TABLE_PI + "{population}/{chromosome}.snps"),
        vs = temp(TABLE_PI + "{population}/{chromosome}.tsv"),
    params:
        mincount = config["popoolation_params"]["tajimad"]["mincount"],
        mincoverage = config["popoolation_params"]["tajimad"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["tajimad"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["tajimad"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["tajimad"]["poolsize"],
        stepsize = config["popoolation_params"]["tajimad"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimad"]["windowsize"],
    log: TABLE_PI + "{population}/{chromosome}.log"
    benchmark: TABLE_PI + "{population}/{chromosome}.json"
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


rule tajimapi_merge_vs:
    input:
        expand(
            TABLE_PI + "{population}/{chromosome}.tsv",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_PI + "{population}.tsv.gz")
    log: PLOT_PI + "merge_vs.log"
    benchmark: PLOT_PI + "merge_vs.json"
    threads: 8
    shell:
        "(bash src/variance_sliding_to_genomic_score.sh {input} "
        "| pigz --best --processes {threads} > {output}) "
        "2> {log}"


rule tajimapi_merge_snps:
    input:
        expand(
            TABLE_PI + "{population}/{chromosome}.snps",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_PI + "{population}.snps.gz")
    threads: 8
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"



rule tajimapi_plot_population:
    """
    Plot a genome-wide Tajima's D distribution
    """
    input:
        tsv_gz = PLOT_PI + "{population}.tsv.gz"
    output:
        z_pdf = protected(PLOT_PI + "{population}_z.pdf"),
        pdf = protected(PLOT_PI + "{population}.pdf")
    params:
        merged_tsv = PLOT_PI + "{population}.tsv"
    log: PLOT_PI + "{population}.log"
    benchmark: PLOT_PI + "{population}.json"
    shell:
        "Rscript src/plot_score.R "
            "--input {input.tsv_gz} "
            "--output {output.pdf} "
        "2>> {log} ; "
        "Rscript src/plot_score.R "
            "--normalize "
            "--input {input.tsv_gz} "
            "--output {output.z_pdf} "
        "2>> {log}"
