include: "__functions__.smk"
include: "mpileup.smk"
include: "variance_sliding.smk"
include: "plot.smk"

rule popoolation:
    input:
        rules.popoolation__mpileup.input,
        rules.popoolation__variance_sliding.input,
        rules.popoolation__plot.input,
