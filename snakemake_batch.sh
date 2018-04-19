#!/usr/bin/env bash
set -euxo pipefail

snakemake \
    --jobs 999 \
    --cluster-config cluster.json \
    --cluster "sbatch --cpus-per-task {threads} --mem {cluster.mem} --time {cluster.time}"
