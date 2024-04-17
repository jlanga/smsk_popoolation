rule popoolation2__plot__merge__:
    input:
        tsvs=lambda w: [
            POP2_FST / f"{chromosome}.w{w.window}-s{w.step}.tsv.gz"
            for chromosome in CHROMOSOMES
        ],
    output:
        tsv_gz=POP2_PLOTS / "all.w{window}-s{step}.tsv.gz",
    log:
        POP2_PLOTS / "all.w{window}-s{step}.log",
    conda:
        "__environment__.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule popoolation2__plot__split_pair__:
    """Split fst table into a pair of populations"""
    input:
        merged_tsv_gz=POP2_PLOTS / "all.w{window}-s{step}.tsv.gz",
    output:
        fst_tsv=POP2_PLOTS / "{pop1}_{pop2}.w{window}-s{step}.fst.gz",
    log:
        POP2_PLOTS / "{pop1}_{pop2}.w{window}-s{step}.fst.log",
    params:
        pop1_number=lambda w: POPULATIONS.index(w.pop1) + 1,
        pop2_number=lambda w: POPULATIONS.index(w.pop2) + 1,
    conda:
        "__environment__.yml"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.merged_tsv_gz} \
        | python3 workflow/scripts/fst_to_genomic_score.py \
            {params.pop1_number} \
            {params.pop2_number} \
        | gzip \
        > {output.fst_tsv} \
        ) 2> {log}
        """


rule popoolation2__plot__:
    """Plot pairwise F_ST distributions over a genome"""
    input:
        fst_tsv=POP2_PLOTS / "{pop1}_{pop2}.w{window}-s{step}.fst.gz",
    output:
        pdf=POP2_PLOTS / "{pop1}_{pop2}.w{window}-s{step}.pdf",
    log:
        POP2_PLOTS / "{pop1}_{pop2}.w{window}-s{step}.pdf.log",
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


rule popoolation2__plot:
    """Make every plot"""
    input:
        [
            POP2_PLOTS / f"{pop1}_{pop2}.w{window}-s{step}.pdf"
            for pop1 in POPULATIONS
            for pop2 in POPULATIONS[POPULATIONS.index(pop1) + 1 :]
            for window, step in POP2_WINDOW_STEP
        ],
