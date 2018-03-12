rule tajimapi_popoolation:
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
        poolsize = lambda wildcards: config["pool_sizes"][wildcards.population],
        stepsize = config["popoolation_params"]["tajimad"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimad"]["windowsize"],
    log: TABLE_PI + "{population}/{chromosome}.log"
    benchmark: TABLE_PI + "{population}/{chromosome}.json"
    conda: "tajimapi.yml"
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
    Merge per chromosome Variance Sliding results into a single file.
    """
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
    conda: "tajimapi.yml"
    shell:
        "(bash src/variance_sliding_to_genomic_score.sh {input} "
        "| pigz --best --processes {threads} > {output}) "
        "2> {log}"


rule tajimapi_merge_snps:
    """
    Merge per chromosome results into a single file.
    """
    input:
        expand(
            TABLE_PI + "{population}/{chromosome}.snps",
            chromosome = CHROMOSOMES,
            population = ["{population}"]
        )
    output: protected(PLOT_PI + "{population}.snps.gz")
    threads: 8
    conda: "tajimapi.yml"
    shell: "pigz --best --keep --stdout --processes {threads} {input} > {output}"



rule tajimapi_plot:
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
    conda: "tajimapi.yml"
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
