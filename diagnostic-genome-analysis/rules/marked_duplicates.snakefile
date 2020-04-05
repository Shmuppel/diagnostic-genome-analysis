rule mark_duplicates:
    """
    Marks duplicates that may have arisen by PCR artifacts.
    
    Input: 
        A sorted BAM file containing the mapped samples.
    
    Output: 
        A new BAM file, in which duplicates have been identified.
        A metrics file indicating the numbers of duplicates.
    
    Shell clarification: 
        java -jar <path to picard.jar> 
        Markduplicates I=<input file> O=<output file path> M=<metrics file>
        ASSUME_SORTED=<true/false> (true as BAM is sorted by samtools sort)
    
    Citation: 
        “Picard Toolkit.” 2019. Broad Institute, GitHub Repository.
        http://broadinstitute.github.io/picard/; Broad Institute
    """
    input:
        tool = config["picard_tool"],
        sample = "runs/{sample}/temp_files/samtools_sort_{sample}.bam"
    output:
        bam = temp("runs/{sample}/temp_files/marked_duplicates_{sample}.bam"),
        metrics = temp("runs/{sample}/temp_files/marked_duplicates_{sample}_metrics.txt")
    benchmark:
        "runs/{sample}/benchmarks/marked_duplicates.txt"
    shell:
        "java -jar {input.tool} "
        "MarkDuplicates I={input.sample} O={output.bam} M={output.metrics} "
        "ASSUME_SORTED=true"