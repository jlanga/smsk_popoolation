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
