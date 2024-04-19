# map
def compose_rg_tag(wildcards):
    identifier = f"ID:{wildcards.population}_{wildcards.library}"
    library = f"LB:truseq_{wildcards.library}"
    platform = "PL:Illumina"
    sample = f"SM:{wildcards.population}"
    return f"@RG\\t{identifier}\\t{library}\\t{platform}\\t{sample}"


## mpileup
def get_library_files_from_sample(wildcards):
    """TODO: needs improvement/simplification
    Return the list of libraries corresponding to a population and chromosome.
    """
    chromosome = wildcards.chromosome
    population = wildcards.population
    libraries = samples[samples["population"] == population]["library"].values.tolist()
    files = [
        PRE_SPLIT / f"{population}.{library}" / f"{chromosome}.bam"
        for library in libraries
    ]
    return files


# coverage
def get_files_for_jellyfish(wildcards):
    population = wildcards.population
    files = [
        READS / f"{population}.{library}_{end}.fq.gz"
        for library in samples[samples["population"] == population]["library"]
        for end in ["1", "2"]
    ]
    return files
