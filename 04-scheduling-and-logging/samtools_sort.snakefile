rule samtools_sort:
    message: "Sorting reads > Executing samtools sort on the following {input} to generate the following {output}."
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input}"
        " > {output}"