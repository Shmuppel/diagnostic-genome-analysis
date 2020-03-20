rule mark_duplicates:
    """
    Marks duplicates that may have arisen by PCR artifacts.
    
    Input: A sorted BAM file containing the mapped samples.
    
    Output: A new BAM file, in which duplicates have been identified.
        A metrics file indicating the numbers of duplicates.
    
    Shell: java -jar <path to picard.jar> 
        Markduplicates I=<input BAM file> O=<output BAM file> M=<metrics file>
        ASSUME_SORTED=true (as BAM is sorted by samtools sort)
    
    Citation: “Picard Toolkit.” 2019. Broad Institute, GitHub Repository.
        http://broadinstitute.github.io/picard/; Broad Institute
    """
    input:
        tool = config["picard_tool"],
        sample = "sorted_reads/{sample}.bam"
    output:
        bam = "marked_duplicates/marked_{sample}.bam",
        metrics = "marked_duplicates/marked_{sample}_metrics.txt"
    shell:
        "java -jar {input.tool} "
        "MarkDuplicates I={input.sample} O={output.bam} M={output.metrics} "
        "ASSUME_SORTED=true"