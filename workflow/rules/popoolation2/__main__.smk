include: "__functions__.smk"
include: "sync.smk"
include: "fst_sliding.smk"


rule popoolation2:
    input:
        rules.popoolation2__sync.input,
        rules.popoolation2__fst_sliding.input,
