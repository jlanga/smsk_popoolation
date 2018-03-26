rule theta_popoolation:
    """
    Get the sliding Tajima's D values.
    """
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{population}.{chromosome}.mpileup.gz"
    output:
        snps = temp(TABLE_THETA + "{population}/{population}.{chromosome}.snps"),
        vs = temp(TABLE_THETA + "{population}/{population}.{chromosome}.tsv"),
    params:
        mincount = config["popoolation_params"]["tajimad"]["mincount"],
        mincoverage = config["popoolation_params"]["tajimad"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["tajimad"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["tajimad"]["mincoveredfraction"],
        poolsize = lambda wildcards: config["pool_sizes"][wildcards.population],
        stepsize = config["popoolation_params"]["tajimad"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimad"]["windowsize"],
    log: TABLE_THETA + "{population}/{population}.{chromosome}.log"
    benchmark: TABLE_THETA + "{population}/{population}.{chromosome}.json"
    conda: "theta.yml"
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
    Merge theta results and take the genomic score
    """
    input:
        expand(
            TABLE_THETA + "{population}/{population}.{chromosome}.tsv",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_THETA + "{population}.tsv.gz")
    log: PLOT_THETA + "merge_vs.log"
    benchmark: PLOT_THETA + "merge_vs.json"
    threads: 8
    conda: "theta.yml"
    shell:
        "(bash src/variance_sliding_to_genomic_score.sh {input} "
        "| pigz --best --processes {threads} > {output}) "
        "2> {log}"


rule theta_merge_snps:
    input:
        expand(
            TABLE_THETA + "{population}/{population}.{chromosome}.snps",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_THETA + "{population}.snps.gz")
    threads: 8
    conda: "theta.yml"
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"



rule theta_plot:
    """
    Plot a genome-wide Tajima's D distribution
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
    conda: "theta.yml"
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
