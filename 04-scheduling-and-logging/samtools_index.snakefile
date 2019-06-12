rule samtools_index:
    message: "Indexing sorted aligned reads > Executing samtools index on the following {input} to generate the following {output}."
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"