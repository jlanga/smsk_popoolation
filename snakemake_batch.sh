#!/usr/bin/env bash
set -euxo pipefail

snakemake \
    --jobs 999 \
    --use-conda \
    --cluster-config cluster.json \
    --cluster "sbatch --job-name {rule} --cpus-per-task {threads} --mem {cluster.mem} --time {cluster.time}"
