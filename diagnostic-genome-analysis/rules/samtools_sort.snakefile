rule samtools_sort:
    """
    Sorts aligned reads by leftmost coordinates.
    
    Input: 
        BAM file to sort.
    
    Output: 
        Sorted BAM file.
        
    Threads:
        6, this number of threats provided a good maximum beyond which paralellization did 
        not significantly improve performance. Depending on sample sizes used it might be
        effective to increase or decrease the number of threads.
        
    Shell clarification: 
        samtools sort -O <output format> <input file> > <output file path>
    
    Citation: 
        Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., 
        Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009)
        The Sequence alignment/map (SAM) format and SAMtools. 
        Bioinformatics, 25, 2078-9. [PMID: 19505943]
    """
    input:
        "runs/{sample}/temp_files/bwa_mem_{sample}.bam.fq.gz"
    output:
        temp("runs/{sample}/temp_files/samtools_sort_{sample}.bam")
    benchmark:
        "runs/{sample}/benchmarks/samtools_sort.txt"
    log:
        "runs/{sample}/logs/samtools_sort.log"
    threads: min(6, workflow.cores - config["reserve_annovar_db_thread"])
    shell:
        "(samtools sort -@ {threads} -O bam {input} > {output}) "
        "2> {log}"