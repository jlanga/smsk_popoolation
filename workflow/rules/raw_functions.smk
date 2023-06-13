def get_reads(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    forward_, reverse_ = samples[
        (samples["population"] == pop) & (samples["library"] == lib)
    ][["forward", "reverse"]].values.tolist()[0]
    return forward_, reverse_
