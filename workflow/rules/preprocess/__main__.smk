include: "__functions__.smk"
include: "map.smk"
include: "mpileup.smk"


rule preprocess:
    input:
        rules.preprocess__mpileup.input,
