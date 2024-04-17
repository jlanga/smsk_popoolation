rule popoolation2__fst_sliding__compute__:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a true file
    """
    input:
        sync=POP2_SUB / "{chromosome}.sync",
    output:
        tsv_gz=POP2_FST / "{chromosome}.w{window}-s{step}.tsv.gz",
    params:
        tsv=lambda w: POP2_FST / f"{w.chromosome}.w{w.window}-s{w.step}.tsv",
        window_size=lambda w: w.window,
        step_size=lambda w: w.step,
        min_covered_fraction=get_min_covered_fraction,
        min_coverage=get_min_coverage,
        max_coverage=get_max_coverage,
        pool_size=compose_population_sizes,
        min_count=get_min_count,
    log:
        POP2_FST / "{chromosome}.w{window}-s{step}.log",
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


rule popoolation2__fst_sliding:
    input:
        [
            POP2_FST / f"{chromosome}.w{window}-s{step}.tsv.gz"
            for chromosome in CHROMOSOMES
            for window, step in POP2_WINDOW_STEP
        ],
