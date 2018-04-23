#!/usr/bin/env bash
set -euxo pipefail

snakemake \
    --jobs 999 \
    --use-conda \
    --cluster-config cluster.json \
    --cluster "sbatch --job-name {rule} --ntasks {threads} --mem {cluster.mem} --time {cluster.time}"
