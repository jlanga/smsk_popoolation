# Popoolation
def get_popoolation_min_count(wildcards):
    return params["popoolation"][wildcards.analysis]["min_count"]


def get_popoolation_min_covered_fraction(wildcards):
    return params["popoolation"][wildcards.analysis]["min_covered_fraction"]


def get_popoolation_step_size(wildcards):
    return params["popoolation"][wildcards.analysis]["step_size"]


def get_popoolation_window_size(wildcards):
    return params["popoolation"][wildcards.analysis]["window_size"]


def get_popoolation_min_coverage(wildcards):
    return params["popoolation"][wildcards.analysis]["min_coverage"]


def get_pool_size(wildcards):
    return (
        samples[samples["population"] == wildcards.population][["pool_size"]]
        .drop_duplicates()
        .values.tolist()[0][0]
    )


def get_popoolation_max_coverage(wildcards):
    return (
        samples[samples["population"] == wildcards.population][["max_coverage"]]
        .drop_duplicates()
        .values.tolist()[0]
    )


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
