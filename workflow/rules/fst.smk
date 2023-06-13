def get_window_size(wildcards):
    return params["popoolation2"]["fst"]["window_size"]


def get_step_size(wildcards):
    return params["popoolation2"]["fst"]["step_size"]


def get_min_covered_fraction(wildcards):
    return params["popoolation2"]["fst"]["min_covered_fraction"]


def get_min_coverage(wildcards):
    return params["popoolation2"]["fst"]["min_coverage"]


def get_max_coverage(wildcards):
    return samples["max_coverage"].max()


def get_min_count(wildcards):
    return params["popoolation2"]["fst"]["min_count"]


def compose_population_sizes(wildcards):
    pool_sizes = (
        samples
        [["population", "pool_size"]]
        .drop_duplicates()
        ["pool_size"]
        .values
        .tolist()
    )
    return ":".join(map(str, pool_sizes))


rule fst_sliding:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a true file
    """
    input:
        sync = SYNC_SUB + "{chromosome}.sync"
    output:
        tsv = temp(TABLE_FST + "{chromosome}.tsv"),
        tsv_gz = protected(TABLE_FST + "{chromosome}.tsv.gz")
    params:
        sync = SYNC_SUB + "{chromosome}.sync",
        tsv = TABLE_FST + "{chromosome}.tsv",
        window_size = get_window_size,
        step_size = get_step_size,
        min_covered_fraction = get_min_covered_fraction,
        min_coverage = get_min_coverage,
        max_coverage = get_max_coverage,
        pool_size = compose_population_sizes,
        min_count = get_min_count
    log:
        TABLE_FST + "{chromosome}.log"
    benchmark:
        TABLE_FST + "{chromosome}.json"
    conda:
        "../envs/fst.yml"
    shell:
        """
        perl workflow/scripts/popoolation2_1201/fst-sliding.pl \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --suppress-noninformative \
            --input {params.sync} \
            --min-covered-fraction {params.min_covered_fraction} \
            --min-coverage {params.min_coverage} \
            --max-coverage {params.max_coverage} \
            --min-count {params.min_count} \
            --output {params.tsv} \
            --pool-size {params.pool_size} \
        2> {log} 1>&2
        gzip --best --keep {params.tsv} 2>> {log}
        """

rule fst_merge:
    input:
        tsvs = expand(
            TABLE_FST + "{chromosome}.tsv",
            chromosome=CHROMOSOMES
        )
    output:
        tsv_gz = protected(
            PLOT_FST + "all.tsv.gz"
        )
    log:
        PLOT_FST + "merge.log"
    benchmark:
        PLOT_FST + "merge.json"
    threads:
        24
    conda:
        "../envs/fst.yml"
    shell:
        "pigz --best --stdout {input} > {output}"


rule fst_split_table:
    """Split fst table into a pair"""
    input:
        merged_tsv_gz = PLOT_FST + "all.tsv.gz"
    output:
        fst_tsv = PLOT_FST + "{pop1}_{pop2}.fst.tsv"
    # log:
    #     PLOT_FST + "split_{pop1}_{pop2}.log"
    # benchmark:
    #     PLOT_FST + "split_{pop1}_{pop2}.json"
    params:
        pop1 = "{pop1}",
        pop2 = "{pop2}"
    conda:
        "../envs/fst.yml"
    shell:
        """
        (gzip --decompress --stdout {input.merged_tsv_gz} \
        | python3 workflow/scripts/fst_to_genomic_score.py \
            {params.pop1} \
            {params.pop2} \
        > {output.fst_tsv})
        """


rule fst_plot:
    """Plot pairwise F_ST distributions over a genome"""
    input:
        fst_tsv = PLOT_FST + "{pop1}_{pop2}.fst.tsv"
    output:
        pdf = PLOT_FST + "{pop1}_{pop2}.pdf"
    log:
        PLOT_FST + "plot_{pop1}_{pop2}.log"
    benchmark:
        PLOT_FST + "plot_{pop1}_{pop2}.json"
    conda:
        "../envs/fst.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input.fst_tsv} \
            --output {output.pdf} \
            --normalize \
        2> {log} 1>&2
        """


rule fst:
    """Make every plot"""
    input:
        [
            PLOT_FST + str(i) + "_" + str(j) + ".pdf"
            for i in range(1, len(POPULATIONS))
            for j in range(i + 1, len(POPULATIONS) + 1)
        ]
