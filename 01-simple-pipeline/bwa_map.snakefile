rule bwa_map:
    message: "Mapping samples against reference genome > Executing bwa mem on the following {input} to generate the following {output}."
    input:
        config["genome"],
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"