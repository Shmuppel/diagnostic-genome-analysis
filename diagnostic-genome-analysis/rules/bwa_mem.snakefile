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
        genome = config["genome"],
        samples = expand("trimmed_samples/{{sample}}_{num}_P.fq.gz",
                         num=["R1","R2"]),
        index = expand("{genome}.{ext}",
                       genome=config["genome"],
                       ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        "mapped_reads/mapped_{sample}.bam"
    threads: 8
    shell:
        "bwa mem -M -t {threads} {input.genome} {input.samples} | samtools view -b > {output}"