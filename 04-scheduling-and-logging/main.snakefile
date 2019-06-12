workdir: "/Users/Shmuppel/Bio-Informatics/data-processing-snakemake/01-simple-pipeline"
configfile: "config.yaml"
include: "bwa_map.snakefile"
include: "samtpols_sort.snakefile"
include: "samtools_index.snakefile"
include: "bcftools_call.snakefile"
include: "report.snakefile"

rule all:
    message: "Executing all snakefile rules to construct a report on the variances of the samples/genome in the 'data' directory."
    input:
        "report.html"