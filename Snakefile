shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


# Folder variables
include: "bin/snakefiles/folders"

# Other variables
POPULATIONS  = config["samples_pe"]
CHROMOSOMES  = config["chromosomes"].split()
ENDS = "1 2 u".split()


# Wildcards functions
def priority_by_chromosome(wildcards):
    priority = 50 - int(wildcards.chromosome) + 1
    return str(priority) 



def tajimad_files(wildcards):
    files = [tajimad + "/" + wildcards.population + "/" + chromosome + ".tsv.gz"
                for chromosome in CHROMOSOMES]
    return files



def tajimapi_files(wildcards):
    files = [tajimapi + "/" + wildcards.population + "/" + chromosome + ".tsv.gz"
                for chromosome in CHROMOSOMES]
    return files



def theta_files(wildcards):
    files = [theta + "/" + wildcards.population + "/" + chromosome + ".tsv.gz"
                for chromosome in CHROMOSOMES]
    return files



rule all:
    input:
        RAW + "genome.fa",
        expand(
            QC + "{population}_{end}.fq.gz",
            population = POPULATIONS,
            end = ENDS
        ),
        expand(
            MAP_FILT + "{population}/{chromosome}.bam",
            population = POPULATIONS,
            chromosome = CHROMOSOMES
        ),
        expand(
            MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz",
            population = POPULATIONS,
            chromosome = CHROMOSOMES
        ),
        expand(
            PLOT_D + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_PI + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_T + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_HP + "{population}.pdf",
            population = POPULATIONS
        ),
        [ 
            PLOT_FST + str(i) + "_" + str(j) +".pdf" 
            for i in range(1, len(POPULATIONS))
            for j in range(i+1, len(POPULATIONS)+1)
        ]
        


rule clean:
    """
    rm -rf data/fastq_raw
    rm -rf data/fastq_trimmed
    rm -rf data/mapping
    """


include: "bin/snakefiles/raw"
include: "bin/snakefiles/qc"
include: "bin/snakefiles/map"
include: "bin/snakefiles/mpileup"
include: "bin/snakefiles/tajimad"
include: "bin/snakefiles/tajimapi"
include: "bin/snakefiles/theta"
include: "bin/snakefiles/hp"
include: "bin/snakefiles/sync"
include: "bin/snakefiles/fst"






    
