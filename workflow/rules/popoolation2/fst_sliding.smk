rule popoolation2__fst_sliding__compute__:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a true file
    """
    input:
        sync=POP2_SUB / "{chromosome}.sync",
    output:
        tsv_gz=POP2_FST / "{chromosome}.tsv.gz",
    params:
        tsv=lambda w: POP2_FST / f"{w.chromosome}.tsv",
        window_size=get_window_size,
        step_size=get_step_size,
        min_covered_fraction=get_min_covered_fraction,
        min_coverage=get_min_coverage,
        max_coverage=get_max_coverage,
        pool_size=compose_population_sizes,
        min_count=get_min_count,
    log:
        POP2_FST / "{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation2/fst-sliding.pl \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --suppress-noninformative \
            --input {input.sync} \
            --min-covered-fraction {params.min_covered_fraction} \
            --min-coverage {params.min_coverage} \
            --max-coverage {params.max_coverage} \
            --min-count {params.min_count} \
            --output {params.tsv} \
            --pool-size {params.pool_size} \
        2> {log} 1>&2

        gzip {params.tsv} 2>> {log}
        """


rule popoolation2__fst_sliding__merge__:
    input:
        tsvs=[POP2_FST / f"{chromosome}.tsv.gz" for chromosome in CHROMOSOMES],
    output:
        tsv_gz=POP2_PLOTS / "all.tsv.gz",
    log:
        POP2_PLOTS / "merge.log",
    threads: 24
    conda:
        "__environment__.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule popoolation2__fst_sliding:
    input:
        POP2_PLOTS / "all.tsv.gz",
