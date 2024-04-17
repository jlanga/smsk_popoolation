rule popoolation__variance_sliding__:
    """
    Run popoolation's Variance sliding script: Tajima's D, Tajima's Theta or Pi
    """
    input:
        mpileup_gz=POP1_SUB / "{population}.{chromosome}.mpileup.gz",
    output:
        snps=temp(POP1_VS / "{population}.{chromosome}.{analysis}.snps"),
        snps_gz=POP1_VS / "{population}.{chromosome}.{analysis}.snps.gz",
        vs=temp(POP1_VS / "{population}.{chromosome}.{analysis}.tsv"),
        vs_gz=POP1_VS / "{population}.{chromosome}.{analysis}.tsv.gz",
    params:
        measure=lambda w: w.analysis,
        min_count=get_popoolation_min_count,
        min_coverage=get_popoolation_min_coverage,
        max_coverage=get_popoolation_max_coverage,
        min_covered_fraction=get_popoolation_min_covered_fraction,
        pool_size=get_pool_size,
        step_size=get_popoolation_step_size,
        window_size=get_popoolation_window_size,
    log:
        POP1_VS / "{population}.{chromosome}.{analysis}.log",
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
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --input <(gzip --decompress --stdout {input.mpileup_gz}) \
            --output {output.vs} \
            --snp-output {output.snps} \
        2> {log} 1>&2

        gzip --keep {output.snps} {output.vs} 2>> {log} 1>&2
        """


rule popoolation__variance_sliding:
    input:
        [
            POP1_VS / f"{population}.{chromosome}.{analysis}.{extension}.gz"
            for population in POPULATIONS
            for chromosome in CHROMOSOMES
            for analysis in ["D", "pi", "theta"]
            for extension in ["snps", "tsv"]
        ],
