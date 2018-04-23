rule popoolation_variance_sliding:
    """Run popoolation's Variance sliding script: Tajima's D, Tajima's Theta or Theta"""
    input:
        mpileup_gz = MPILEUP_SUB + "{population}/{population}.{chromosome}.mpileup.gz"
    output:
        snps = temp(TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.snps"),
        vs = temp(TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.tsv"),
    params:
        measure = "{analysis}",
        mincount = lambda wildcards: config["popoolation_params"][wildcards.analysis]["min_count"],
        mincoverage = lambda wildcards: config["popoolation_params"][wildcards.analysis]["min_coverage"],
        maxcoverage = lambda wildcards: config["samples"][wildcards.population]["max_coverage"],
        mincoveredfraction = lambda wildcards: config["popoolation_params"][wildcards.analysis]["min_covered_fraction"],
        poolsize = lambda wildcards: config["samples"][wildcards.population]["pool_size"],
        stepsize = lambda wildcards: config["popoolation_params"][wildcards.analysis]["step_size"],
        windowsize = lambda wildcards: config["popoolation_params"][wildcards.analysis]["window_size"],
    log: TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.log"
    benchmark: TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.json"
    conda: "popoolation.yml"
    shell:
        "perl src/popoolation_1.2.2/Variance-sliding.pl "
            "--measure {params.measure} "
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



rule popoolation_merge_variance_sliding:
    """Merge results across chromosomes"""
    input:
        expand(
            TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.tsv",
            chromosome = CHROMOSOMES,
            population = ["{population}"],
            analysis = ["{analysis}"]
        )
    output: protected(PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.tsv.gz")
    log: PLOT_POPOOLATION + "merge_vs.log"
    benchmark: PLOT_POPOOLATION + "merge_vs.json"
    threads: 8
    conda: "popoolation.yml"
    shell:
        "(bash src/variance_sliding_to_genomic_score.sh {input} "
        "| pigz --best --processes {threads} > {output}) "
        "2> {log}"


rule popoolation_merge_snps:
    """Merge snps"""
    input:
        expand(
            TABLE_POPOOLATION + "{analysis}/{population}.{chromosome}.{analysis}.snps",
            chromosome = CHROMOSOMES,
            population = ["{population}"],
            analysis = ["{analysis}"]
        )
    output: protected(PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.snps.gz")
    threads: 8
    conda: "popoolation.yml"
    shell: "pigz --processes {threads} --best --stdout --processes {threads} {input} > {output}"



rule popoolation_plot:
    """Plot a genome-wide result's from variance sliding"""
    input:
        tsv_gz = PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.tsv.gz"
    output:
        pdf = protected(PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.pdf")
    log: PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.log"
    benchmark: PLOT_POPOOLATION + "{analysis}/{population}.{analysis}.json"
    conda: "popoolation.yml"
    shell:
        "Rscript src/plot_score.R "
            "--input {input.tsv_gz} "
            "--output {output.pdf} "
        "2> {log}"



rule popoolation_d:
    input:
        expand(
            PLOT_POPOOLATION + "D/{population}.D.pdf",
            population= POPULATIONS
        )


rule popoolation_pi:
    input:
        expand(
            PLOT_POPOOLATION + "pi/{population}.pi.pdf",
            population= POPULATIONS
        )


rule popoolation_theta:
    input:
        expand(
            PLOT_POPOOLATION + "theta/{population}.theta.pdf",
            population= POPULATIONS
        )
