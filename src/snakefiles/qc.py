def get_adaptor(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    adaptor = (
        samples
        [(samples["population"] == pop) & (samples["library"] == lib)]
        ["adaptor"]
        .values
        .tolist()
        [0]
    )
    return adaptor


def get_phred(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    phred = (
        samples
        [(samples["population"] == pop) & (samples["library"] == lib)]
        ["phred"]
        .values
        .tolist()
        [0]
    )
    return phred


def get_trimmomatic_params(wildcards):
    # The YAML file introduces a \n
    return params["trimmomatic"]["extra"].strip()


rule qc_trimmomatic:
    """Run trimmomatic on paired end mode

    to eliminate Illumina adaptors andremove low quality regions and reads.

    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from
    Trimmomatic.

    Further on they are catted and compressed with gzip/pigz (level 1).
    Sequences will be stored permanently later on on CRAM
    """
    input:
        fwd = RAW + "{population}.{library}_1.fq.gz",
        rev = RAW + "{population}.{library}_2.fq.gz"
    output:
        fwd = temp(QC + "{population}.{library}_1.fq.gz"),
        rev = temp(QC + "{population}.{library}_2.fq.gz"),
        fwd_unp = temp(QC + "{population}.{library}_3.fq.gz"),
        rev_unp = temp(QC + "{population}.{library}_4.fq.gz")
    params:
        adaptor = get_adaptor,
        phred = get_phred,
        trimmomatic_params = get_trimmomatic_params
    log:
        QC + "{population}.{library}.trimmomatic_pe.log"
    benchmark:
        QC + "{population}.{library}.trimmomatic_pe.json"
    threads:
        4
    priority:
        50  # Do this and later the mappings
    conda:
        "qc.yml"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -{params.phred} \
            <(gzip --decompress --stdout {input.fwd}) \
            <(gzip --decompress --stdout {input.rev}) \
            >(gzip --fast > {output.fwd}) \
            >(gzip --fast > {output.fwd_unp}) \
            >(gzip --fast > {output.rev}) \
            >(gzip --fast > {output.rev_unp}) \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log} 1>&2
        """

rule qc:
    input:
        [
            QC + population + "." + library + "_1.fq.gz"
            for population, library in (
                samples
                [["population", "library"]]
                .values
                .tolist()
            )
        ]
