rule bwa_mem:
    """
    Performs local alignment of samples to a reference genome.

    Input: 
        Indexed reference genome files and trimmed paired end reads. 

    Output: 
        A BAM file containing the mapped sample.
        
    Threads: 
        12
    
    Shell clarification: 
        bwa mem -M (for picard compatibility) {input genome} {input sample reads} | 
        samtools view -b (to export into BAM file) {output file path}

    Reference & further info: 
        http://bio-bwa.sourceforge.net/bwa.shtml

    Citation: 
        Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. 
        arXiv:1303.3997v1 [q-bio.GN].
    """
    input:
        genome = config["genome"],
        index = expand("{genome}.{ext}",
                       genome=config["genome"],
                       ext=["amb", "ann", "bwt", "pac", "sa"]),
        samples = expand("runs/{{sample}}/temp_files/trimmomatic_{read}_P.fq.gz",
                         read=["R1", "R2"])
    output:
        temp("runs/{sample}/temp_files/bwa_mem_{sample}.bam.fq.gz")
    benchmark:
        "runs/{sample}/benchmarks/bwa_mem.txt"
    threads: min(12, workflow.cores - config["reserve_annovar_db_thread"])
    shell:
        "bwa mem -M -t {threads} {input.genome} {input.samples} | samtools view -b > {output}"