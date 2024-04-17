rule popoolation2__mpileup__join__:
    """
    Join multiple single-sample mpileups into one
    """
    input:
        mpileups=lambda w: [
            PRE_MPILEUP / f"{population}.{w.chromosome}.mpileup.gz"
            for population in POPULATIONS
        ],
    output:
        temp(POP2_MPILEUP / "{chromosome}.mpileup.gz"),
    log:
        POP2_MPILEUP / "{chromosome}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( python workflow/scripts/join_mpileups.py \
            {input.mpileups} \
        | gzip \
        > {output} \
        ) 2> {log}
        """


rule popoolation2__mpileup:
    input:
        [POP2_MPILEUP / f"{chromosome}.mpileup.gz" for chromosome in CHROMOSOMES],
