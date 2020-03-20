rule annovar:
	"""
	Input: vcf file containing filtered variants.

	Output:

	Shell: -remove <remove temp files>
        -nastring <fill fields without annotiation>

	Citation: Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic
		variants from next-generation sequencing data, Nucleic Acids Research,
		38:e164, 2010
	"""
	input:
		tool=config["annovar_tool"],
		humandb=config["humandb"],
		vcf="vcf/vcf_A.vcf"
	output:
		"annovar/{sample}.hg19_multianno.txt",
		"annovar/{sample}.hg19_multianno.vcf"
	run:
		# Download the appropiate databases
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene {input.humandb} ")
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 {input.humandb} ")
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp138 {input.humandb} ")
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all {input.humandb} ")
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct {input.humandb} ")
		shell("{input.tool}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20190305 {input.humandb} ")
		# Perform table.annover.pl
		shell("{input.tool}table_annovar.pl -vcfinput {input.vcf} -out annovar/{wildcards.sample} -buildver hg19 {input.humandb} "
			  "-polish -remove "
			  "-protocol refGene,exac03,avsnp138,ljb26_all,ALL.sites.2014_10,clinvar_20190305 "
			  "-operation g,f,f,f,f,f ")
		# Clean up databases
		# shell("rm {input.humandb}*")