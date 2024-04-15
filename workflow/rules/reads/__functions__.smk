def get_forward(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    forward_ = samples[(samples["population"] == pop) & (samples["library"] == lib)][
        ["forward"]
    ].values.tolist()[0]
    return forward_


def get_reverse(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    reverse_ = samples[(samples["population"] == pop) & (samples["library"] == lib)][
        ["reverse"]
    ].values.tolist()[0]
    return reverse_
