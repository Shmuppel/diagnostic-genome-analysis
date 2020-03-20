def get_sample_from_wildcard(wildcards):
    """
    Retrieves the files associated from a sample stored in config.yaml.
    Which sample/file to retrieve is based on the wildcard passed as an argument to this function
    either manually or dynamically through 'input = get_sample_from_wildcard'.

    This function (as does the whole pipeline) assumes 2 paired end reads files per sample are stored
    with file name formatting: '{sample}_R1.fastq', ''{sample}_R2.fastq'.
    """
    return expand(
        f"{config['samples'][wildcards.sample]}_{{num}}.fastq",
        num=["R1","R2"])

rule fastqc:
    """
    Perform FastQC Analysis on a sample.

    Input: 
        Single End FASTQ sample(s).

    Output: 
        Analysis of samples, summarized in an HTML file. 
        An archieve of all graphs and figures.

    Shell clarification: 
        cat {input file} | fastqc --outdir=<output file directory> stdin:{output file path}

    Reference & further info:
        https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

    Citation: 
        Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
        Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
    """
    input:
        sample = get_sample_from_wildcard
    output:
        "fastqc/{sample}_fastqc.html"
    shell:
        "cat {input} | fastqc --outdir=./fastqc stdin:{wildcards.sample}"
