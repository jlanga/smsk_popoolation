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
        | awk
            '{{hist[$4]++}} END {{for (i in hist) print "{params.population}" "\\t" i "\\t" hist[i]}}' \
        | sort -n -k2,2 \
        > {output.hist} \
        ) 2> {log}
        """


rule preprocess__coverage:
    input:
        [PRE_COV / f"{population}.hist" for population in POPULATIONS],
