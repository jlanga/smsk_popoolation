shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"


singularity: "docker://continuumio/miniconda3:4.4.10"

# Folder variables
include: "src/snakefiles/folders.py"

# Other variables
POPULATIONS_PE = config["samples_pe"] if config["samples_pe"] is not None else []
POPULATIONS_SE = config["samples_se"] if config["samples_se"] is not None else []
POPULATIONS = [x for x in POPULATIONS_PE] + [x for x in POPULATIONS_SE]
PAIRS = ["pe_pe", "pe_se"]
CHROMOSOMES  = config["chromosomes"].split(" ")
ENDS = "1 2 u".split(" ")


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

include: "src/snakefiles/generic.py"
include: "src/snakefiles/raw.py"
include: "src/snakefiles/qc.py"
include: "src/snakefiles/map.py"
include: "src/snakefiles/mpileup.py"
include: "src/snakefiles/tajimad.py"
include: "src/snakefiles/tajimapi.py"
include: "src/snakefiles/theta.py"
include: "src/snakefiles/hp.py"
include: "src/snakefiles/sync.py"
include: "src/snakefiles/fst.py"
include: "src/snakefiles/reports.py"


rule all:
    input:
        rules.qc.input,
        rules.map.input,
        # raw rules
        # RAW + "genome.fa",
        # trimming
        # expand(
        #    QC + "{population}.{library}_{end}.fq.gz",
        #    population = POPULATIONS,
        #    end = "1 2".split()
        # ),
        # mapping
        #expand(
        #    MAP_FILT + "{population}/{chromosome}.bam",
        #    population = POPULATIONS,
        #    chromosome = CHROMOSOMES
        #),
        #expand(
        #    MPILEUP_SUB + "{population}/{chromosome}.mpileup.gz",
        #    population = POPULATIONS,
        #    chromosome = CHROMOSOMES
        #),
        expand(
            PLOT_D + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_PI + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_THETA + "{population}.pdf",
            population = POPULATIONS
        ),
        expand(
            PLOT_HP + "{population}.pdf",
            population = POPULATIONS
        ),
        rules.fst.input,
        rules.reports.input



rule clean:
    shell:
        "rm -r results/"
