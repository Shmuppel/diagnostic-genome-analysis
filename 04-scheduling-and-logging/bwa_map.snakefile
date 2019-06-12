rule bwa_map:
    message: "Mapping samples against reference genome > Executing bwa mem on the following {input} to generate the following {output}."
    threads: 8
    input:
        config["genome"],
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"