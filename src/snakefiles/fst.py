rule fst_sliding_chromosome:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a
    """
    input:
        sync = SYNC_SUB + "{chromosome}.sync"
    output:
        tsv = temp(
            TABLE_FST + "{chromosome}.tsv"
        ),
        tsv_gz = protected(
            TABLE_FST + "{chromosome}.tsv.gz"  # TEMP!
        )
    params:
        sync = SYNC_SUB + "{chromosome}.sync",
        tsv  = TABLE_FST + "{chromosome}.tsv",
        window_size = config["popoolation2_params"]["fst_sliding"]["window_size"],
        step_size = config["popoolation2_params"]["fst_sliding"]["step_size"],
        min_covered_fraction = config["popoolation2_params"]["fst_sliding"]["min_covered_fraction"],
        min_coverage = config["popoolation2_params"]["fst_sliding"]["min_coverage"],
        max_coverage = 10000,
        # max_coverage = ",".join([
        #     config["samples"][population]["max_coverage"]
        #     for population in POPULATIONS
        # ]),
        pool_size = ":".join([
            config["samples"][population]["pool_size"]
            for population in POPULATIONS
        ]),
        min_count = config["popoolation2_params"]["fst_sliding"]["min_count"]
    log: TABLE_FST + "{chromosome}.log"
    benchmark: TABLE_FST + "{chromosome}.json"
    conda: "fst.yml"
    shell:
        "perl src/popoolation2_1201/fst-sliding.pl "
            "--window-size {params.window_size} "
            "--step-size {params.step_size} "
            "--suppress-noninformative "
            "--input {params.sync} "
            "--min-covered-fraction {params.min_covered_fraction} "
            "--min-coverage {params.min_coverage} "
            "--max-coverage {params.max_coverage} "
            "--min-count {params.min_count} "
            "--output {params.tsv} "
            "--pool-size {params.pool_size} "
        "2> {log} 1>&2 ; "
        "gzip --best --keep {params.tsv} 2>> {log}"



rule fst_merge:
    input:
        tsvs= expand(
            TABLE_FST + "{chromosome}.tsv",
            chromosome = CHROMOSOMES
        )
    output:
        tsv_gz = protected(
            PLOT_FST + "all.tsv.gz"
        )
    log: TABLE_FST + "merge.log"
    benchmark: TABLE_FST + "merge.json"
    threads: 24
    conda: "fst.yml"
    shell:
        "pigz --best --stdout {input} > {output}"


rule fst_split_table:
    """Split fst table into a pair"""
    input:
        merged_tsv_gz = PLOT_FST + "all.tsv.gz"
    output:
        fst_tsv = PLOT_FST + "{pop1}_{pop2}.fst.tsv"
    log:
        PLOT_FST + "split_{pop1}_{pop2}.log"
    benchmark:
        PLOT_FST + "split_{pop1}_{pop2}.json"
    params:
        pop1 = "{pop1}",
        pop2 = "{pop2}"
    threads:
        1
    conda: "fst.yml"
    shell:
        "(gzip --decompress --stdout {input.merged_tsv_gz} "
        "| python3 src/fst_to_genomic_score.py "
            "{params.pop1} "
            "{params.pop2} "
        "> {output.fst_tsv}) "
        "2> {log} 1>&2"


rule fst_plot:
    """Plot pairwise F_ST distributions over a genome"""
    input:
        fst_tsv = PLOT_FST + "{pop1}_{pop2}.fst.tsv"
    output:
        pdf = PLOT_FST + "{pop1}_{pop2}.pdf"
    threads: 1
    log: PLOT_FST + "plot_{pop1}_{pop2}.log"
    benchmark: PLOT_FST + "plot_{pop1}_{pop2}.json"
    conda: "fst.yml"
    shell:
        "Rscript src/plot_score.R "
            "--input {input.fst_tsv} "
            "--output {output.pdf} "
        "2> {log} 1>&2; "


rule fst:
    """Make every plot"""
    input:
        [
            PLOT_FST + str(i) + "_" + str(j) +".pdf"
            for i in range(1, len(POPULATIONS))
            for j in range(i+1, len(POPULATIONS)+1)
        ]
