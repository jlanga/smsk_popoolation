def get_window_size(wildcards):
    return params["popoolation2"]["fst"]["window_size"]


def get_step_size(wildcards):
    return params["popoolation2"]["fst"]["step_size"]


def get_min_covered_fraction(wildcards):
    return params["popoolation2"]["fst"]["min_covered_fraction"]


def get_min_coverage(wildcards):
    return params["popoolation2"]["fst"]["min_coverage"]


def get_max_coverage(wildcards):
    return samples["max_coverage"].max()


def get_min_count(wildcards):
    return params["popoolation2"]["fst"]["min_count"]


def compose_population_sizes(wildcards):
    pool_sizes = (
        samples[["population", "pool_size"]]
        .drop_duplicates()["pool_size"]
        .values.tolist()
    )
    return ":".join(map(str, pool_sizes))
