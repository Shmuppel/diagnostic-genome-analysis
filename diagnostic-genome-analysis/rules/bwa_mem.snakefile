rule bwa_mem:
    """
    Performs local alignment of samples to a reference genome.

    Input: 
        Indexed reference genome. 
        Trimmed paired end reads. 

    Output: 
        A BAM file containing the mapped sample.

    Shell clarification: 
        bwa mem -M (for picard compatibility) <input genome> <input sample reads> | 
        samtools view -b (to export into BAM file) <output file path>

    Reference & further info: 
        http://bio-bwa.sourceforge.net/bwa.shtml

    Citation: 
        Li H. (2013) Aligning sequence reads, clone sequences and
        assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN].
    """
    input:
        config["genome"],
        "trimmed_samples/{sample}_R1_P.fq.gz",
        "trimmed_samples/{sample}_R2_P.fq.gz"
    output:
        "mapped_reads/mapped_{sample}.bam"
    threads: 8
    shell:
        "bwa mem -M -t {threads} {input} | samtools view -b > {output}"