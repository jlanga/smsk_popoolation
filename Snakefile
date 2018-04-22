shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"

singularity: "docker://continuumio/miniconda3:4.4.10"

# Folder variables
include: "src/snakefiles/folders.py"

# Other variables
POPULATIONS = config["samples"] if config["samples"] is not None else []
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
# include: "src/snakefiles/tajimad.py"
# include: "src/snakefiles/tajimapi.py"
# include: "src/snakefiles/theta.py"
include: "src/snakefiles/popoolation.py"
include: "src/snakefiles/hp.py"
include: "src/snakefiles/sync.py"
include: "src/snakefiles/fst.py"
include: "src/snakefiles/reports.py"


rule all:
    input:
        # rules.qc.input,
        # rules.map.input,
        # rules.mpileup_subsampled.input,
        rules.popoolation_d.input,
        rules.popoolation_pi.input,
        rules.popoolation_theta.input,
        expand(
            PLOT_HP + "{population}.pdf",
            population = POPULATIONS
        ),
        rules.fst.input,
        rules.reports.input



rule clean:
    shell:
        "rm -r results/"
