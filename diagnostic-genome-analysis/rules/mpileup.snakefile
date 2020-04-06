rule samtools_mpileup:
    """
    Convert BAM file into a pileup format file.
    
    Input: 
        Reference genome and BAM file to convert.
    
    Output: 
        File in mpileup format containing reads.
        
    Threads:
        1, as mpileup does not support multithreading.
    
    Shell clarification: 
        samtools mpileup -f <reference genome> <input file> -o <output file path>
    
    Citation: 
        Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., 
        Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009)
        The Sequence alignment/map (SAM) format and SAMtools. 
        Bioinformatics, 25, 2078-9. [PMID: 19505943]
    """
    input:
        reference_genome = config["genome"],
        bam = "runs/{sample}/temp_files/split_reads/{sample}.REF_chr{chr}.bam"
    output:
        mpileup = temp("runs/{sample}/temp_files/mpileup_chr{chr}.mpileup")
    threads: 1
    benchmark:
        "runs/{sample}/benchmarks/samtools_mpileup_{chr}.txt"
    log:
        "runs/{sample}/logs/mpileup_{chr}.log"
    shell:
        "(samtools mpileup -f {input.reference_genome} {input.bam} > {output.mpileup}) "
        "2> {log}"