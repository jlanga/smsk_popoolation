rule popoolation2__fst_sliding__compute__:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a true file
    """
    input:
        sync=POP2_SUB / "{chromosome}.sync",
    output:
        tsv=temp(POP2_TABLES / "{chromosome}.tsv"),
        tsv_gz=POP2_TABLES / "{chromosome}.tsv.gz",
    params:
        sync=POP2_SUB / "{chromosome}.sync",
        tsv=POP2_TABLES / "{chromosome}.tsv",
        window_size=get_window_size,
        step_size=get_step_size,
        min_covered_fraction=get_min_covered_fraction,
        min_coverage=get_min_coverage,
        max_coverage=get_max_coverage,
        pool_size=compose_population_sizes,
        min_count=get_min_count,
    log:
        POP2_TABLES / "{chromosome}.log",
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
            --output {output.tsv} \
            --pool-size {params.pool_size} \
        2> {log} 1>&2

        gzip --keep {output.tsv} 2>> {log}
        """


rule popoolation2__fst_sliding__merge__:
    input:
        tsvs=expand(POP2_TABLES / "{chromosome}.tsv", chromosome=CHROMOSOMES),
    output:
        tsv_gz=POP2_PLOTS / "all.tsv.gz",
    log:
        POP2_PLOTS / "merge.log",
    threads: 24
    conda:
        "__environment__.yml"
    shell:
        "pigz --stdout {input} > {output}"


rule popoolation2__fst_sliding__split__:
    """Split fst table into a pair"""
    input:
        merged_tsv_gz=POP2_PLOTS / "all.tsv.gz",
    output:
        fst_tsv=POP2_PLOTS / "{pop1}_{pop2}.fst.tsv.gz",
    log:
        POP2_PLOTS / "split_{pop1}_{pop2}.log",
    params:
        pop1="{pop1}",
        pop2="{pop2}",
    conda:
        "__environment__.yml"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.merged_tsv_gz} \
        | python3 workflow/scripts/fst_to_genomic_score.py \
            {params.pop1} \
            {params.pop2} \
        | gzip \
        > {output.fst_tsv} \
        ) 2> {log}
        """


rule popoolation2__fst_sliding__plot__:
    """Plot pairwise F_ST distributions over a genome"""
    input:
        fst_tsv=POP2_PLOTS / "{pop1}_{pop2}.fst.tsv.gz",
    output:
        pdf=POP2_PLOTS / "{pop1}_{pop2}.pdf",
    log:
        POP2_PLOTS / "plot_{pop1}_{pop2}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_score.R \
            --input {input.fst_tsv} \
            --output {output.pdf} \
            --normalize \
        2> {log} 1>&2
        """


rule popoolation2__fst_sliding:
    """Make every plot"""
    input:
        [
            POP2_PLOTS / f"{str(i)}_{str(j)}.pdf"
            for i in range(1, len(POPULATIONS))
            for j in range(i + 1, len(POPULATIONS) + 1)
        ],
