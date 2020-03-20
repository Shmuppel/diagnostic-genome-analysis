def get_sample_from_wildcard(wildcards):
    """
    Retrieves the files associated with a sample stored in config.yaml.

    Which sample/file to retrieve is based on the wildcard passed as an argument
    to this function dynamically through 'input = get_sample_from_wildcard'.

    Assumes 2 paired end reads files per sample are stored with file name formatting:
    '{sample}_R1.fastq', ''{sample}_R2.fastq'.
    """
    return expand(
        f"{config['samples'][wildcards.sample]}_{{num}}.fastq",
        num=["R1","R2"])

rule trimmomatic:
    """ 
    Trims reads, removes ILLUMINACLIP artifacts, and removes low length reads.

    Input: 
        Two FASTAQ files. Trimmomatic is set to use paired end trimming
        both forward and reverse reads are expected to be present.

    Output: 
        Trimmed sequence reads. Reads that did not make the cut have been 
        moved to unpaired archieves (1U and 2U).
    
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
        sample = get_sample_from_wildcard
    output:
        expand(["trimmed_samples/{{sample}}_R1_{replicate}.fq.gz",
                "trimmed_samples/{{sample}}_R2_{replicate}.fq.gz"],
                replicate=["P", "U"])
    shell:
        "java -jar {input.tool} "
        "PE {input.sample} {output} "
        "ILLUMINACLIP:tools/trimmomatic-0.32-1/share/TruSeq3-PE.fa:2:30:10 "
        "SLIDINGWINDOW:4:20 "
        "MINLEN:70"
