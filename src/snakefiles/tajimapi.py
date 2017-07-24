rule tajimapi_table_population_chromosome:
    """
    Get the sliding Tajima's Pi values.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps = temp(TABLE_PI + "{population}/{chromosome}.snps"),
        vs = temp(TABLE_PI + "{population}/{chromosome}.tsv"),
    params:
        mincount = config["popoolation_params"]["tajimapi"]["mincount"],
        mincoverage = config["popoolation_params"]["tajimapi"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["tajimapi"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["tajimapi"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["tajimapi"]["poolsize"],
        stepsize = config["popoolation_params"]["tajimapi"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimapi"]["windowsize"],
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
    """
    Merge all tsvs into a single tsv.gz
    """
    input:
        expand(
            TABLE_PI + "{population}/{chromosome}.tsv",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_PI + "{population}.tsv.gz")
    threads: 8
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"


rule tajimapi_merge_snps:
    """
    Merge all snps into a single snps.tsv
    """
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
    Plot a genome-wide Tajima's Pi distribution
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
            "none "
            "{input.tsv_gz} "
            "{output.pdf} "
        "2>> {log} ; "
        "Rscript src/plot_score.R "
            "z "
            "{input.tsv_gz} "
            "{output.z_pdf} "
        "2>> {log}"
