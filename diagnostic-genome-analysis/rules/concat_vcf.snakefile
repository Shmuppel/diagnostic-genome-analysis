rule concat_vcf:
    """
    Helper rule that concatenates all VCF files generated for a sample into a single VCF file.
    
    Input:
        Sorted VCF files to be concatenated.
        
    Output:
        Single VCF containing all significant variants of the sample. 
    """
    input:
        expand("runs/{{sample}}/temp_files/vcf_chr{chr}.mpileup", chr=config["chromosomes"])
    output:
        temp("runs/{sample}/temp_files/vcf_{sample}.vcf"),
    benchmark:
        "runs/{sample}/benchmarks/concat_vcf.txt"
    shell:
        "cat {input} > {output}"