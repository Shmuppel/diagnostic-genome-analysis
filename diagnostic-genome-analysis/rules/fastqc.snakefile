def get_sample_from_wildcard(wildcards):
    """
    Retrieves the files associated from a sample stored in config.yaml.
    Which sample/file to retrieve is based on the wildcard passed as an argument to this function
    either manually or dynamically through 'input = get_sample_from_wildcard'.

    This function (as does the whole pipeline) assumes 2 paired end reads files per sample are
    stored with file name formatting: '{sample}_R1.fastq', ''{sample}_R2.fastq'.
    """
    return expand("{sample}_{num}.fastq",
                  num=["R1","R2"],
                  sample=config['samples'][wildcards.sample])

rule fastqc:
    """
    Perform FastQC Analysis on a sample's read files.

    Input: 
        Single end read FASTQ files.

    Output: 
        Analysis of samples, summarized in an HTML file. 
        An archieve of all graphs and figures.
        
    Threads: 
        2, this rule will issue two FastQC processes running in parallel, each analyzing 1 
        single end read file. 
        
        FastQC natively supports multithreading, however that caused some complications
        with Snakemake dependencies due to FastQC's multiple output files and their naming
        convention. The current GNU parallel command isn't much different in practice, and
        allows renaming of the output files into a more accessible naming convention. 

    Shell clarification: 
        parallel -j {no. threads to use}
        cat {input file} | fastqc --outdir=<output file directory> stdin:{sample name}_R<job no.>

    Reference & further info:
        https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

    Citation: 
        Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
        Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
    """
    input:
        sample = get_sample_from_wildcard
    output:
        expand("results/fastqc/{{sample}}_{num}_fastqc.html", num=["R1", "R2"])
    threads: min(2, workflow.cores - config["reserve_annovar_db_thread"])
    shell:
        "parallel -j {threads} " \
        "'cat {{1}} | fastqc --outdir=./results/fastqc stdin:{wildcards.sample}_R{{#}}' ::: {input} "
