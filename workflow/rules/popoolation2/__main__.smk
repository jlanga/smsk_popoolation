include: "__functions__.smk"
include: "mpileup.smk"
include: "sync.smk"
include: "fst_sliding.smk"
include: "plot.smk"


rule popoolation2:
    input:
        rules.popoolation2__sync.input,
        rules.popoolation2__fst_sliding.input,
        rules.popoolation2__plot.input,
