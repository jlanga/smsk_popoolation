include: "__functions__.smk"


rule reads__link__:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{population}.{library}_1.fq.gz",
        reverse_=READS / "{population}.{library}_2.fq.gz",
    log:
        READS / "{population}.{library}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2> {log}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2>> {log}
        """


rule reads:
    input:
        [
            READS / f"{population}.{library}_{end}.fq.gz"
            for population, library in POPULATION_LIBRARY
            for end in [1, 2]
        ],
