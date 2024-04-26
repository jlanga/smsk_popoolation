rule preprocess__coverage__compute_hist__:
    input:
        mpileups=lambda w: [
            PRE_MPILEUP / w.population / f"{chromosome}.mpileup.gz"
            for chromosome in CHROMOSOMES
        ],
    output:
        hist=PRE_COV / "{population}.hist",
    log:
        PRE_COV / "{population}.hist.log",
    conda:
        "__environment__.yml"
    params:
        population=lambda w: w.population,
    shell:
        """
        ( gzip -dc {input.mpileups} \
        | awk \
            '{{hist[$4]++}} END {{for (i in hist) print "{params.population}" "\\t" i "\\t" hist[i]}}' \
        | sort -n -k2,2 \
        > {output.hist} \
        ) 2> {log}
        """


rule preprocess_coverage__plot__:
    input:
        [PRE_COV / f"{population}.hist" for population in POPULATIONS],
    output:
        pdf=PRE_COV / "coverage.pdf",
        coverage=PRE_COV / "coverage.tsv",
    log:
        PRE_COV / "coverage.log",
    conda:
        "__environment__.yml"
    shell:
        """
        Rscript workflow/scripts/plot_coverage.R 2> {log} 1>&2
        """


rule preprocess__coverage:
    input:
        hist=PRE_COV / "coverage.pdf",
        coverage=PRE_COV / "coverage.tsv",
