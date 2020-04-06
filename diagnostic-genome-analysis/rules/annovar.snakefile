rule annovar:
	"""
	Annotates genetic variants using various databases, which can be configured in config.yaml.
	
	Input:
		Annovar annotation databases to perform operations with, retrieved from get_remote_annovar_dbs.
		VCF file containing variants to annotate. 
		
	Output:
		Multiannotated .txt and .vcf files containing annotated variants.
		
	Params:
		out_dir: path pointing to the output directory.
		protocol: databases to annotate with, retrieved from config.yaml.
		operation: operation to perform for each database, retrieved from config.yaml.
	
	Shadow: 
        Full, creates a temporary environment to store the output values of Annovar.
        Facilitates automatic removal of unwanted / unused files. 
	
	Threads:
		4, Annovar support multithreading to some extend. More than 4 threads showed
		no significant improvements in runtime. 
		
	Shell clarification:
		<annovar tool path>table_annovar.pl 
		-thread <number of threads to use>
		-vcfinput <input VCF file> 
		-out <output file location>
		-buildver <genome build version>
		-protocol <databases to annotate with>
		-operation <operation to perform for each database> <path to Annovar databases>
		-remove (remove temporary files when done)
		
	Citation: 
		Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic
		variants from next-generation sequencing data, Nucleic Acids Research,
		38:e164, 2010
	"""
	input:
		rules.get_remote_annovar_dbs.output,
		vcf = "runs/{sample}/temp_files/vcf_{sample}.vcf"
	output:
		# wildcard sample used dynamically with config.yaml; genome build is added as a prefix by Annovar.
		expand(["runs/{{sample}}/results/{{sample}}.{genome_build}_multianno.txt",
				"runs/{{sample}}/results/{{sample}}.{genome_build}_multianno.vcf"],
		   		genome_build = config["genome_build"])
	params:
		out_dir = "runs/{sample}/results/{sample}",
		# gather string of databases to use for variant annotating.
		protocols = ",".join(config["annovar_dbs"].keys()),
		# gather string of operations to perform with these databases.
		operations = ",".join([config["annovar_dbs"][db]["operation"] for db in config["annovar_dbs"]])
	shadow: "full"
	benchmark:
		"runs/{sample}/benchmarks/annovar.txt"
	log:
		"runs/{sample}/logs/annovar.log"
	threads: min(4, workflow.cores)
	shell:
		"({config[annovar_tool]}table_annovar.pl "
		"-thread {threads} "
		"-vcfinput {input.vcf} "
		"-out {params.out_dir} "
		"-buildver {config[genome_build]} "
		"-protocol {params.protocols} "
		"-operation {params.operations} "
		"{config[annovar_db_storage]} "
		"-remove ) "
		"2> {log}"
