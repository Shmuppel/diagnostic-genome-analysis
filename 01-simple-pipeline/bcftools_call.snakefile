rule bcftools_call:
    message: "Extracting variants from aggregated mapped reads > Executing mpileup and bcftools on the following {input} to generate the following {output}."
    input:
        fa= config["genome"],
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"