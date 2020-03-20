rule samtools_mpileup:
    """
    Convert BAM file into a pileup format file.
    
    Input: 
        Reference genome and BAM file to convert.
    
    Output: 
        Pileup file containing samples.
    
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
        bam = "marked_duplicates/marked_{sample}.bam",
    output:
        mpileup = "mpileup/mpileup_{sample}.mpileup"
    shell:
        "samtools mpileup -f {input.reference_genome} {input.bam} -o {output.mpileup}"