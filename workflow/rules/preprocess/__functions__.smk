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
        MAP_FILT / f"{population}.{library}.{chromosome}.cram" for library in libraries
    ]
    return files


def get_indel_window(wildcards):
    return params["popoolation"]["find_indels"]["indel_window"]


def get_indel_min_count(wildcards):
    return params["popoolation"]["find_indels"]["min_count"]


def get_subsample_min_qual(wildcards):
    return params["popoolation"]["subsample"]["min_qual"]


def get_subsample_method(wildcards):
    return params["popoolation"]["subsample"]["method"]


def get_subsample_max_coverage(wildcards):
    return samples[samples["population"] == wildcards.population][
        "max_coverage"
    ].tolist()[0]


def get_subsample_target_coverage(wildcards):
    return params["popoolation"]["subsample"]["target_coverage"]
