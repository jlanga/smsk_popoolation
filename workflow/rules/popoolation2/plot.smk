rule popoolation2__plot__split_pair__:
    """Split fst table into a pair of populations"""
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


rule popoolation2__plot__:
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


rule popoolation2__plot:
    """Make every plot"""
    input:
        [
            POP2_PLOTS / f"{str(i)}_{str(j)}.pdf"
            for i in range(1, len(POPULATIONS))
            for j in range(i + 1, len(POPULATIONS) + 1)
        ],
