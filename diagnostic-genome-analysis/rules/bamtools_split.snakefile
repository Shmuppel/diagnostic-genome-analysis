rule bamtools_split:
    """
    Splits the reads of a sample into separate bam files, one file per chromosome.
    
    Input: 
        The mapped, marked reads of a sample contained in a BAM file.
    
    Output: 
        Multiple BAM files, depending on the amount of chromosomes specified in config.yaml.
        Each file contains the reads of one chromosome.
    
    Params:
        out_dir: the path pointing towards the output directory, used as a stub by bamtools split.
        
    Shadow: 
        Full, creates a temporary environment to store the output values of bamtools split.
        Allows unused chromosome seperated BAM files to be discarded.
    
    Shell clarification:
        bamtools split -in <input>
        -reference (to separate by chromosome) 
        -stub <prefix/directory of output path>"
    
    Citation: 
        Derek W. Barnett, Erik K. Garrison, Aaron R. Quinlan, Michael P. Strömberg, Gabor T. Marth, 
        BamTools: a C++ API and toolkit for analyzing and managing BAM files, Bioinformatics, 
        Volume 27, Issue 12, 15 June 2011, Pages 1691–1692, https://doi.org/10.1093/bioinformatics/btr174
    """
    input:
        "runs/{sample}/temp_files/marked_duplicates_{sample}.bam"
    output:
        temp(expand("runs/{{sample}}/temp_files/split_reads/{{sample}}.REF_chr{chr}.bam",
                    chr = config["chromosomes"]))
    params:
        out_dir = "runs/{sample}/temp_files/split_reads/{sample}"
    shadow: "full"
    wildcard_constraints: chr="."
    benchmark:
        "runs/{sample}/benchmarks/bamtools_split.txt"
    log:
        "runs/{sample}/logs/bamtools_split.log"
    shell:
        "(bamtools split -in {input} -reference -stub {params.out_dir}) "
        "2> {log}"