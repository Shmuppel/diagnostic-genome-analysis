rule bwa_mem:
    """
    Performs local alignment of samples to a reference genome.

    Input: Reference genome, indexed; samples that have undergone
        trimming.

    Output: One BAM file containing the mapped sample.

    Shell: bwa mem -M (for picard compatibility) <input genome and samples> |
        samtools view -b (to export into BAM file) <output path>

    Reference & further info: 
        http://bio-bwa.sourceforge.net/bwa.shtml

    Citation: Li H. (2013) Aligning sequence reads, clone sequences and
        assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN].
    """
    input:
        config["genome"],
        "trimmed_samples/{sample}_R1_P.fq.gz",
        "trimmed_samples/{sample}_R2_P.fq.gz"
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -M -t {threads} {input} | samtools view -b > {output}"