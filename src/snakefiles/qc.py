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
    Number of threads used:
        4 for trimmomatic
        2 for gzip inputs
        2 for gzip outputs
        Total: 8
    """
    input:
        forward = RAW + "{sample}_1.fq.gz",
        reverse = RAW + "{sample}_2.fq.gz"
    output:
        forward     = temp(QC + "{sample}_1.fq.gz"),
        reverse     = temp(QC + "{sample}_2.fq.gz"),
        unpaired    = protected(QC + "{sample}.final.pe_se.fq.gz")
    params:
        unpaired_1  = QC + "{sample}_3.fq.gz",
        unpaired_2  = QC + "{sample}_4.fq.gz",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log:
        QC + "trimmomatic_pe_{sample}.log"
    benchmark:
        QC + "trimmomatic_pe_{sample}.json"
    threads:
        24 # I've been able to work with pigz and 24 trimmomatic threads.
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -{params.phred} \
            <(pigz -dc {input.forward} ) \
            <(pigz -dc {input.reverse} ) \
            >(cut -f 1 -d " " | pigz -9 > {output.forward} ) \
            {params.unpaired_1} \
            >(cut -f 1 -d " " | pigz -9 > {output.reverse} ) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}

        zcat {params.unpaired_1} {params.unpaired_2} |
        cut -f 1 -d " " |
        pigz -9 > {output.unpaired}

        rm {params.unpaired_1} {params.unpaired_2}
        """



rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and
        remove low quality regions and reads.
    Input is piped through gzip/pigz.
    Output is piped to gzip.
    Threads used:
        4 for trimmomatic
        1 for gzip input
        1 for gzip output
    """
    input:
        single = RAW + "{sample}_se.fq.gz",
    output:
        single = protected(QC + "{sample}.final.se.fq.gz")
    params:
        adaptor = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    log:
        QC + "trimmomatic_se_{sample}.log"
    benchmark:
        QC + "trimmomatic_se_{sample}.json"
    threads:
        8
    shell:
        """
        trimmomatic SE \
            -threads {threads} \
            -{params.phred} \
            <(pigz -dc {input.single}) \
            >(cut -f 1 -d " " | pigz -9 > {output.single}) \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}
        """