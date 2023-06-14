rule fst_sliding:
    """
    Compute sliding F_STs in one chromosome.

    Note: fst-sliding requires a true file
    """
    input:
        sync=SYNC_SUBSAMPLED / "{chromosome}.sync",
    output:
        tsv=temp(FST_TABLES / "{chromosome}.tsv"),
        tsv_gz=protected(FST_TABLES / "{chromosome}.tsv.gz"),
    params:
        sync=SYNC_SUBSAMPLED / "{chromosome}.sync",
        tsv=FST_TABLES / "{chromosome}.tsv",
        window_size=get_window_size,
        step_size=get_step_size,
        min_covered_fraction=get_min_covered_fraction,
        min_coverage=get_min_coverage,
        max_coverage=get_max_coverage,
        pool_size=compose_population_sizes,
        min_count=get_min_count,
    log:
        FST_TABLES / "{chromosome}.log",
    conda:
        "../envs/fst.yml"
    shell:
        """
        perl workflow/scripts/popoolation2_1201/fst-sliding.pl \
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
        gzip --best --keep {output.tsv} 2>> {log}
        """


rule fst_merge:
    input:
        tsvs=expand(FST_TABLES / "{chromosome}.tsv", chromosome=CHROMOSOMES),
    output:
        tsv_gz=protected(FST_PLOTS / "all.tsv.gz"),
    log:
        FST_PLOTS / "merge.log",
    threads: 24
    conda:
        "../envs/fst.yml"
    shell:
        "pigz --best --stdout {input} > {output}"


rule fst_split_table:
    """Split fst table into a pair"""
    input:
        merged_tsv_gz=FST_PLOTS / "all.tsv.gz",
    output:
        fst_tsv=FST_PLOTS / "{pop1}_{pop2}.fst.tsv",
    log:
        FST_PLOTS / "split_{pop1}_{pop2}.log",
    params:
        pop1="{pop1}",
        pop2="{pop2}",
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
        fst_tsv=FST_PLOTS / "{pop1}_{pop2}.fst.tsv",
    output:
        pdf=FST_PLOTS / "{pop1}_{pop2}.pdf",
    log:
        FST_PLOTS / "plot_{pop1}_{pop2}.log",
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
            FST_PLOTS / f"{str(i)}_{str(j)}.pdf"
            for i in range(1, len(POPULATIONS))
            for j in range(i + 1, len(POPULATIONS) + 1)
        ],
