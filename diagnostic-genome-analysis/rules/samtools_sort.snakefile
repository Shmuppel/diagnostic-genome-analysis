rule samtools_sort:
    """
    Sorts alignments by leftmost coordinates.
    
    Input: BAM file to sort.
    
    Output: Sorted BAM file.
    
    Shell: samtools sort -O bam (output format) <input BAM> 
    
    Citation: Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N.,
        Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009)
        The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
    """
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -O bam {input} > {output}"