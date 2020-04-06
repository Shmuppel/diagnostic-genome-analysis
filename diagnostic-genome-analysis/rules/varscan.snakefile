rule varscan:
	"""
	Identifies statistically signifant SNPs of a mpileup file using VarScan 2.

	Input: 
		File in mpileup format to scan for SNPs.

	Output: 
		VCF file containing filtered variants.
	
	Threads:
		1, as varscan does not support multithreading.
		
	Shell clarification: 
		java -jar <path to varscan jar>
		mpileup2snp <mpileup input file>
		--min-var-freq <minimum allele frequency %>
		--p-value <p-value threshold>
		--output-vcf <output to VCF 0/1>

	Citation: 
		Koboldt, D. C. and Zhang, Q. and Larson, D. E. and Shen, D. and McLellan, M. D. and Lin, L. 
		and Miller, C. A. and Mardis, E. R. and Ding, L. and Wilson, R. K. (2012).
		VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing.
		In Genome Research, 22 (3), pp. 568–576.
	"""
	input:
		tool = config["varscan_tool"],
		mpileup = "runs/{sample}/temp_files/mpileup_chr{chr}.mpileup"
	output:
		vcf = temp("runs/{sample}/temp_files/vcf_chr{chr}.mpileup")
	benchmark:
		"runs/{sample}/benchmarks/varscan_{chr}.txt"
	threads: 1
	shell:
		"java -jar {input.tool} "
		"mpileup2snp {input.mpileup} "
		"--min-var-freq 0.3 "
		"--p-value 0.05 "
		"--output-vcf 1 > {output.vcf}"