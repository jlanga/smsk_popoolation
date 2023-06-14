def get_sync_indel_window(wildcards):
    return params["popoolation2"]["find_indels"]["indel_window"]


def get_sync_min_count(wildcards):
    return params["popoolation2"]["find_indels"]["min_count"]


def compose_mpileups_comma(wildcards):
    chromosome = wildcards.chromosome
    mpileups = [
        str(MPILEUP_RAW / f"{population}/{population}.{chromosome}.mpileup.gz")
        for population in POPULATIONS
    ]
    composed = "{" + ",".join(mpileups) + "}"
    return composed


def get_sync_min_qual(wildcards):
    return params["popoolation2"]["subsample"]["min_qual"]


def get_sync_target_coverage(wildcards):
    return params["popoolation2"]["subsample"]["target_coverage"]


def compose_max_coverages(wildcards):
    coverages = (
        samples[["population", "max_coverage"]]
        .drop_duplicates()["max_coverage"]
        .values.tolist()
    )
    coverages = map(str, coverages)
    return ",".join(coverages)


def get_sync_subsample_method(wildcards):
    return params["popoolation2"]["subsample"]["method"]
