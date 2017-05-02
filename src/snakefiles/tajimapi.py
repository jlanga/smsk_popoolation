rule tajimapi_table_population_chromosome:
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz"
    output:
        snps_gz = TABLE_PI + "{population}/{chromosome}.snps.gz",
        vs_gz = TABLE_PI + "{population}/{chromosome}.tsv.gz"
    params:
        snps = TABLE_PI + "{population}/{chromosome}.snps",
        vs = TABLE_PI + "{population}/{chromosome}.tsv",
        mincount = config["popoolation_params"]["tajimapi"]["mincount"],
        mincoverage = config["popoolation_params"]["tajimapi"]["mincoverage"],
        maxcoverage = config["popoolation_params"]["tajimapi"]["maxcoverage"],
        mincoveredfraction = config["popoolation_params"]["tajimapi"]["mincoveredfraction"],
        poolsize = config["popoolation_params"]["tajimapi"]["poolsize"],
        stepsize = config["popoolation_params"]["tajimapi"]["stepsize"],
        windowsize = config["popoolation_params"]["tajimapi"]["windowsize"],
    threads:
        1
    log:
        TABLE_PI + "{population}/{chromosome}.log"
    benchmark:
        TABLE_PI + "{population}/{chromosome}.json"
    shell:
        "perl src/popoolation_1.2.2/Variance-sliding.pl "
            "--measure pi "
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
        


rule tajimapi_plot_population:
    input:
        tsvs = expand(
            TABLE_PI + "{population}/{chromosome}.tsv.gz",
            population = "{population}",
            chromosome = CHROMOSOMES
        )
    output:
        merged_tsv_gz = PLOT_PI + "{population}.tsv.gz", 
        z_pdf = PLOT_PI + "{population}_z.pdf",
        pdf = PLOT_PI + "{population}.pdf"
    params:
        merged_tsv = PLOT_PI + "{population}.tsv"
    threads:
        1
    log:
        PLOT_PI + "{population}.log"
    benchmark:
        PLOT_PI + "{population}.json"
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
