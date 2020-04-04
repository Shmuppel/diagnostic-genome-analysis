def get_sample_from_wildcard(wildcards):
    """
    Retrieves the files associated from a sample stored in config.yaml.
    Which sample/file to retrieve is based on the wildcard passed as an argument to this function
    either manually or dynamically through 'input = get_sample_from_wildcard'.

    This function (as does the whole pipeline) assumes 2 paired end reads files per sample are stored
    with file name formatting: '{sample}_R1.fastq', ''{sample}_R2.fastq'.
    """
    return expand("{sample}_{num}.fastq",
                  num=["R1","R2"],
                  sample=config['samples'][wildcards.sample])

rule trimmomatic:
    """  
    Trims reads, removes ILLUMINACLIP artifacts, and removes low length reads.

    Input: 
        Two FASTAQ files. Trimmomatic is set to use paired end trimming both forward and reverse reads 
        are expected to be present.
        
    Output: 
        Trimmed sequence reads. Reads that did not make the cut have been moved to unpaired archieves 
        (1U and 2U).
        
    Threads:
        2, although trimmomatic supports multithreading its implementation is not significantly 
        effective past 2 threads. For multiple samples 2 threads per sample revolves to multiple 
        trimmomatic rules being able to run at a relatively fast runtime. When dealing with 
        notably large files, such that runtime exacerbates, splitting the input files into batches
        would be an effective solution.
        
    Shell clarification:
        java -jar <path to trimmomatic.jar>
        PE (paired end) <input file> <output file path>
        
        ILLUMINACLIP (removes Illumina adapters)
            <fasta with illumina adapters>:
            <seed mismatches>:
            <palindrome clip treshold>:
            <simple clip treshold>

        SLIDINGWINDOW (performs a sliding window trimming)
            <number of bases to average across>:
            <average quality required>

        MINLEN (removes reads that fall below a given length)
            <minimum length of reads to be kept>

    Reference & further info: 
        http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

    Citation:
        Bolger, A. M., Lohse, M., & Usadel, B. (2014). 
        Trimmomatic: A flexible trimmer for Illumina Sequence Data.
    """
    input:
        tool = config["trimmomatic_tool"],
        illumina_clip = config["illumina_clip"],
        sample = get_sample_from_wildcard
    output:
        temp(expand("trimmed_samples/{{sample}}_{read}_{replicate}.fq.gz",
                    read=["R1", "R2"],
                    replicate=["P", "U"]))
    threads: min(2, workflow.cores - config["reserve_annovar_db_thread"])
    shell:
        "java -jar {input.tool} "
        "PE -threads {threads} {input.sample} {output} "
        "ILLUMINACLIP:{input.illumina_clip}:2:30:10 "
        "SLIDINGWINDOW:4:20 "
        "MINLEN:70"
