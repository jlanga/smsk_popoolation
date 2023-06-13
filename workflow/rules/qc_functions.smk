def get_adaptor(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    adaptor = samples[(samples["population"] == pop) & (samples["library"] == lib)][
        "adaptor"
    ].values.tolist()[0]
    return adaptor


def get_phred(wildcards):
    pop = wildcards.population
    lib = wildcards.library
    phred = samples[(samples["population"] == pop) & (samples["library"] == lib)][
        "phred"
    ].values.tolist()[0]
    return phred


def get_trimmomatic_params(wildcards):
    # The YAML file introduces a \n
    return params["trimmomatic"]["extra"].strip()
