include: "__functions__.smk"
include: "index.smk"
include: "map.smk"
include: "mpileup.smk"


rule preprocess:
    input:
        rules.preprocess__mpileup.input,
