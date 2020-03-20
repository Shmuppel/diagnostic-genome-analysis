rule varscan:
	"""
	Identifies SNPs of a mpileup file using VarScan 2.

	Input: mpileup file to scan for variants.

	Output: vcf file containing filtered variants.

	Shell: java -jar tools/varscan-v2.3.9/VarScan.jar
		   mpileup2snp <mpileup file>
		   --min-var-freq <minimum allelle frequency>
		   --p-value <p-value treshold>
		   --output-vcf <output to vcf 0/1>

	Citation: Koboldt, D. C. and Zhang, Q. and Larson, D. E. and Shen, D. and McLellan,
		M. D. and Lin, L. and Miller, C. A. and Mardis, E. R. and Ding, L. and Wilson, R. K. (2012).
		VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing.
		In Genome Research, 22 (3), pp. 568–576.
	"""
	input:
		tool=config["varscan_tool"],
		mpileup="mpileup/mpileup_{sample}.mpileup"
	output:
		vcf="vcf/vcf_{sample}.vcf"
	shell:
		"java -jar {input.tool} "
		"mpileup2snp {input.mpileup} --min-var-freq 0.3 --p-value 0.99 "
		"--output-vcf 1 > {output.vcf}"