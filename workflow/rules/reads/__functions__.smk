def get_file(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    files = samples[(samples["population"] == pop) & (samples["library"] == lib)][
        ["forward", "reverse"]
    ].values.tolist()[0]
    return files


def get_forward(wildcards):
    return get_file(wildcards)[0]


def get_reverse(wildcards):
    return get_file(wildcards)[1]
