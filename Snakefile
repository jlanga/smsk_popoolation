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

# Folder variables
include: "src/snakefiles/folders.py"

# Other variables
POPULATIONS = samples["population"].drop_duplicates().values.tolist()
PAIRS = ["pe_pe", "pe_se"]
CHROMOSOMES = features["chromosomes"].split(" ")
ENDS = "1 2 u".split(" ")


include: "src/snakefiles/generic.py"
include: "src/snakefiles/raw.py"
include: "src/snakefiles/qc.py"
include: "src/snakefiles/map.py"
include: "src/snakefiles/mpileup.py"
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
