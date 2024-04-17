# Popoolation - find_indels
POP1_INDEL_WINDOW = params["popoolation"]["find_indels"]["indel_window"]
POP1_INDEL_MIN_COUNT = params["popoolation"]["find_indels"]["min_count"]


# Popoolation - subsample
POP1_SUBSAMPLE_MIN_QUAL = params["popoolation"]["subsample"]["min_qual"]
POP1_SUBSAMPLE_METHOD = params["popoolation"]["subsample"]["method"]
POP1_SUBSAMPLE_TARGET_COVERAGE = params["popoolation"]["subsample"]["target_coverage"]


def get_subsample_max_coverage(wildcards):
    return samples[samples["population"] == wildcards.population][
        "max_coverage"
    ].tolist()[0]


# Popoolation - variance_sliding
POP1_VS_MIN_COUNT = params["popoolation"]["variance_sliding"]["min_count"]
POP1_VS_MIN_COVERAGE = params["popoolation"]["variance_sliding"]["min_coverage"]
POP1_VS_MIN_COVERED_FRACTION = params["popoolation"]["variance_sliding"][
    "min_covered_fraction"
]


def get_vs_pool_size(wildcards):
    return (
        samples[samples["population"] == wildcards.population][["pool_size"]]
        .drop_duplicates()
        .values.tolist()[0][0]
    )


def get_vs_max_coverage(wildcards):
    return (
        samples[samples["population"] == wildcards.population][["max_coverage"]]
        .drop_duplicates()
        .values.tolist()[0]
    )
