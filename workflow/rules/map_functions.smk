def compose_rg_tag(wildcards):
    identifier = f"ID:{wildcards.population}_{wildcards.library}"
    library = f"LB:truseq_{wildcards.library}"
    platform = "PL:Illumina"
    sample = f"SM:{wildcards.population}"
    return f"@RG\t{identifier}\t{library}\t{platform}\t{sample}"
