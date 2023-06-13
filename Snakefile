# pylint: disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version

min_version("5.0")

shell.prefix("set -euo pipefail;")

params = yaml.load(open("params.yml", "r"), Loader=yaml.SafeLoader)
features = yaml.load(open("features.yml", "r"), Loader=yaml.SafeLoader)
samples = pd.read_table("samples.tsv")

singularity: "docker://continuumio/miniconda3:4.4.10"


# Other variables
POPULATIONS = samples["population"].drop_duplicates().values.tolist()
PAIRS = ["pe_pe", "pe_se"]
CHROMOSOMES = features["chromosomes"].split(" ")
ENDS = "1 2 u".split(" ")

include: "src/snakefiles/folders.smk"
include: "src/snakefiles/generic.smk"
include: "src/snakefiles/raw.smk"
include: "src/snakefiles/qc.smk"
include: "src/snakefiles/map.smk"
include: "src/snakefiles/mpileup.smk"
include: "src/snakefiles/popoolation.smk"
include: "src/snakefiles/hp.smk"
include: "src/snakefiles/sync.smk"
include: "src/snakefiles/fst.smk"
include: "src/snakefiles/reports.smk"


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
            population=POPULATIONS
        ),
        rules.fst.input,
        rules.reports.input


rule clean:
    shell:
        """
        if [[ -d results ]]; then
            rm -r results/
        fi
        """
