rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from
    Trimmomatic. Further on they are catted and compressed with gzip/pigz
    (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
    header. It is done posterior to the trimming since the output comes
    slower than the input is read.
    """
    input:
        forward = RAW + "{sample}_1.fq.gz",
        reverse = RAW + "{sample}_2.fq.gz"
    output:
        forward     = temp(QC + "{sample}_1.fq.gz"),
        reverse     = temp(QC + "{sample}_2.fq.gz"),
        forward_unp = temp(QC + "{sample}_3.fq.gz"),
        reverse_unp = temp(QC + "{sample}_4.fq.gz")
    params:
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log: QC + "trimmomatic_pe_{sample}.log"
    benchmark: QC + "trimmomatic_pe_{sample}.json"
    threads: 8
    shell:
        "trimmomatic PE "
            "-threads {threads} "
            "-{params.phred} "
            "<(pigz --decompress --stdout {input.forward}) "
            "<(pigz --decompress --stdout {input.reverse}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| pigz --fast > {output.forward}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| pigz --fast > {output.forward_unp}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| pigz --fast > {output.reverse}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| pigz --fast > {output.reverse_unp}) "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} "
        "2> {log}"



rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and
        remove low quality regions and reads.
    Input is piped through pigz.
    Output is piped to pigz.
    """
    input:
        single = RAW + "{sample}_se.fq.gz",
    output:
        single = temp(QC + "{sample}.final.se.fq.gz")
    params:
        adaptor = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log: QC + "trimmomatic_se_{sample}.log"
    benchmark: QC + "trimmomatic_se_{sample}.json"
    threads: 8
    shell:
        "trimmomatic SE "
            "-threads {threads} "
            "-{params.phred} "
            "<(pigz --decompress --stdout {input.single}) "
            ">(cut --fields 1 --delimiter \" \" "
                "| pigz --fast > {output.single}) "
            "ILLUMINACLIP:{params.adaptor}:2:30:10 "
            "{params.trimmomatic_params} "
        "2> {log}"
