include: "__functions__.smk"
include: "index.smk"
include: "map.smk"
include: "mpileup.smk"
include: "coverage.smk"


rule preprocess:
    input:
        rules.preprocess__mpileup.input,
        rules.preprocess__coverage.input,
