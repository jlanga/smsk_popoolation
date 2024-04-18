rule popoolation__variance_sliding__:
    """
    Run popoolation's Variance sliding script: Tajima's D, Tajima's Theta or Pi
    """
    input:
        mpileup_gz=POP1_SUB / "{chromosome}" / "{population}.mpileup.gz",
    output:
        snps_gz=POP1_VS
        / "{chromosome}"
        / "{population}.{analysis}.w{window}-s{step}.snps.gz",
        vs_gz=POP1_VS
        / "{chromosome}"
        / "{population}.{analysis}.w{window}-s{step}.tsv.gz",
    params:
        snps=lambda w: POP1_VS
        / w.chromosome
        / f"{w.population}.{w.analysis}.w{w.window}-s{w.step}.snps",
        vs=lambda w: POP1_VS
        / w.chromosome
        / f"{w.population}.{w.analysis}.w{w.window}-s{w.step}.tsv",
        measure=lambda w: w.analysis,
        min_count=POP1_VS_MIN_COUNT,
        min_coverage=POP1_VS_MIN_COVERAGE,
        max_coverage=get_vs_max_coverage,
        min_covered_fraction=POP1_VS_MIN_COVERED_FRACTION,
        pool_size=get_vs_pool_size,
        step=lambda w: w.step,
        window=lambda w: w.window,
    log:
        POP1_VS / "{chromosome}" / "{population}.{analysis}.w{window}-s{step}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        perl workflow/scripts/popoolation/Variance-sliding.pl \
            --measure {params.measure} \
            --fastq-type sanger \
            --min-count {params.min_count} \
            --min-coverage {params.min_coverage} \
            --max-coverage {params.max_coverage} \
            --min-covered-fraction {params.min_covered_fraction} \
            --pool-size {params.pool_size} \
            --window-size {params.window} \
            --step-size {params.step} \
            --input <(gzip --decompress --stdout {input.mpileup_gz}) \
            --output {params.vs} \
            --snp-output {params.snps} \
        2> {log} 1>&2

        gzip --keep {params.snps} {params.vs} 2>> {log} 1>&2
        """


rule popoolation__variance_sliding:
    input:
        [
            POP1_VS
            / chromosome
            / f"{population}.{analysis}.w{window}-s{step}.{extension}.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
            for analysis in POP1_ANALYSES
            for extension in ["snps", "tsv"]
            for window, step in POP1_WINDOW_STEP
        ],
