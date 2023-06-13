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
        forward_=RAW / "{population}.{library}_1.fq.gz",
        reverse_=RAW / "{population}.{library}_2.fq.gz",
    output:
        forward_=temp(QC / "{population}.{library}_1.fq.gz"),
        reverse_=temp(QC / "{population}.{library}_2.fq.gz"),
        forward_unp=temp(QC / "{population}.{library}_3.fq.gz"),
        reverse_unp=temp(QC / "{population}.{library}_4.fq.gz"),
    params:
        adaptor=get_adaptor,
        phred=get_phred,
        trimmomatic_params=get_trimmomatic_params,
    log:
        QC / "{population}.{library}.trimmomatic_pe.log",
    benchmark:
        QC / "{population}.{library}.trimmomatic_pe.bmk"
    threads: 4
    priority: 50  # Do this and later the mappings
    conda:
        "../envs/qc.yml"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -{params.phred} \
            <(gzip --decompress --stdout {input.forward_}) \
            <(gzip --decompress --stdout {input.reverse_}) \
            >(gzip --fast > {output.forward_}) \
            >(gzip --fast > {output.forward_unp}) \
            >(gzip --fast > {output.reverse_}) \
            >(gzip --fast > {output.reverse_unp}) \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log} 1>&2
        """


rule qc:
    input:
        [
            QC / f"{population}.{library}_{end}.fq.gz"
            for population, library in (
                samples[["population", "library"]].values.tolist()
            )
            for end in ["1", "2"]
        ],
