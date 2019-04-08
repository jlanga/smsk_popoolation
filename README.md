# smsk_popoolation: A Snakemake pipeline for population genomics

[![Build Status](https://travis-ci.org/jlanga/smsk_popoolation.svg?branch=master)](https://travis-ci.org/jlanga/smsk_popoolation)
[![DOI](https://zenodo.org/badge/76841262.svg)](https://zenodo.org/badge/latestdoi/76841262)

## 1. Description

This is a repo that contains installers and snakemake scripts to execute the pipelines described by Kofler et al. in popoolation 1 and 2:

- QC with Trimmomatic

- Mapping with Bowtie2

- BAM wrangling and SNP calling with samtools

- Population measures with popoolation

- Pairwise comparisons between populations with popoolation2

## 2. First steps

Follow the contents of the `.travis.yml` file:

1. Install ([ana](https://www.continuum.io/downloads)|[mini](http://conda.pydata.org/miniconda.html))conda

2. Clone and install the software

    ```sh
    git clone https://github.com/jlanga/smsk_popoolation.git smsk_popoolation
    cd smsk_popoolation
    snakemake --use-conda --create-envs-only
    ```

3. Run the test dataset:

    ```sh
    snakemake --use-conda -
    ```

4. Modify the following files:

    - `features.yaml`: the path to the genome reference, and the names of every chromosome to be processed.

    - `samples.tsv`: paths and library information of each of the samples.

    - `params.yml`: execution parameters of different tools.

    - `cluster.json`: memory limits for some of the steps-

5. Execute the pipeline:

    ```sh
    snakemake --use-conda -j
    ```

## Representation of the pipeline

![smsk_popoolation pipeline](https://cdn.rawgit.com/jlanga/smsk_popoolation/master/rulegraph.svg)

## Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)

- [popoolation](https://sourceforge.net/p/popoolation/wiki/Main/)

- [popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)

- [samtools](http://www.htslib.org/)
