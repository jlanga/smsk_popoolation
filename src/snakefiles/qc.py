rule qc_trimmomatic_pe:
    """Run trimmomatic on paired end mode

    to eliminate Illumina adaptors andremove low quality regions and reads.

    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from Trimmomatic.

    Further on they are catted and compressed with gzip/pigz (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ header. It is done
    posterior to the trimming since the output comes slower than the input is read.
    """
    input:
        forward = RAW + "{population}.{library}_1.fq.gz",
        reverse = RAW + "{population}.{library}_2.fq.gz"
    output:
        forward     = temp(QC + "{population}.{library}_1.fq.gz"),
        reverse     = temp(QC + "{population}.{library}_2.fq.gz"),
        forward_unp = temp(QC + "{population}.{library}_3.fq.gz"),
        reverse_unp = temp(QC + "{population}.{library}_4.fq.gz")
    params:
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.population][wildcards.library]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.population][wildcards.library]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log: QC + "{population}.{library}.trimmomatic_pe.log"
    benchmark: QC + "{population}.{library}.trimmomatic_pe.json"
    threads: 4
    conda: "qc.yml"
    shell:
        "trimmomatic PE "
            "-threads {threads} "
            "-{params.phred} "
            "<(gzip --decompress --stdout {input.forward}) "
            "<(gzip --decompress --stdout {input.reverse}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| gzip --fast > {output.forward}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| gzip --fast > {output.forward_unp}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| gzip --fast > {output.reverse}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| gzip --fast > {output.reverse_unp}) "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} "
        "2> {log}"



rule qc_trimmomatic_se:
    """Run trimmomatic on single end mode
    to eliminate Illumina adaptors and remove low quality regions and reads.
    Input is piped through gzip.
    Output is piped to gzip.
    """
    input:
        single = RAW + "{population}.{library}_se.fq.gz",
    output:
        single = temp(QC + "{population}.{library}_se.fq.gz")
    params:
        adaptor = lambda wildcards: config["samples_se"][wildcards.population][wildcards.library]["adaptor"],
        phred = lambda wildcards: config["samples_se"][wildcards.population][wildcards.library]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log: QC + "{population}.{library}.trimmomatic_se.log"
    benchmark: QC + "{population}.{library}.trimmomatic_se.json"
    threads: 2
    conda: "qc.yml"
    shell:
        "trimmomatic SE "
            "-threads {threads} "
            "-{params.phred} "
            "<(gzip --decompress --stdout {input.single}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| gzip --fast > {output.single}) "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} "
        "2> {log}"


def get_qc_files(wildcards):
    """ TODO: needs improvement/simplification
    Return the list of libraries corresponding to a population and chromosome.
    """
    return [
        QC + population + "." + library + "_1.fq.gz"
        for population in config["samples_pe"]
        for library in config["samples_pe"][population]
    ]


rule qc:
    input:
        get_qc_files
