rule bwa_index:
    """
    Index reference genome in the FASTA format.

    Input: 
        Reference Genome (FASTA Format or compressed FASTA).

    Output: 
        Various files, some in binary, that are used for indexing 
        the reference genome. 
    
    Threads:
        1, BWA index does not support multithreading, therefore it is limited
        to use 1 thread.
        
    Shell: 
        bwa index -a <algorithm for constructing human BWT index> <input file>

    Reference & further info: 
        http://bio-bwa.sourceforge.net/bwa.shtml

    Citation: 
        Li H. (2013) Aligning sequence reads, clone sequences and assembly 
        contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN].
    """
    input:
	    config["genome"]
    output:
        expand("{genome}.{ext}",
               genome=config["genome"],
               ext=["amb", "ann", "bwt", "pac", "sa"])
    threads: 1
    priority: 1
    shell:
        "bwa index -a bwtsw {input}"
