def get_annovar_db_paths():
	"""
    Retrieves the local paths where Annovar databases are/should be stored.

    Used so that databases only have to be downloaded once, and so that dependencies
    can be dynamically added / removed in config.yaml.
    """
	# example db path: annovar databases folder + genome build prefix + annovar database file; for each database.
	return [config["annovar_db_storage"] + config["genome_build"] + '_' + config["annovar_dbs"][db]["file"]
			for db in config["annovar_dbs"]]

rule annovar:
	"""
	Annotates genetic variants using various databases, which can be configured
	in config.yaml.
	
	Input:
		Annovar annotation databases to perform operations with.
		VCF file containing variants to annotate. 
		
	Output:
		Multiannotated .txt and .vcf files containing annotated variants.
		
	Shell clarification:
		<annovar tool path>table_annovar.pl 
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
		get_annovar_db_paths(),
		vcf="vcf/vcf_{sample}.vcf"
	output:
		# wildcard sample used dynamically with config.yaml; genome build is added as a prefix by Annovar.
		expand(["annovar/{{sample}}.{genome_build}_multianno.txt",
				"annovar/{{sample}}.{genome_build}_multianno.vcf"],
		   		genome_build = config["genome_build"])
	run:
		# gather string of databases to use for variant annotating.
		protocols = ",".join(config["annovar_dbs"].keys())
		# gather string of operations to perform with these databases.
		operations = ",".join([config["annovar_dbs"][db]["operation"] for db in config["annovar_dbs"]])

		shell("{config[annovar_tool]}table_annovar.pl "
			  "-vcfinput {input.vcf} " \
			  "-out annovar/{wildcards.sample} "
			  "-buildver {config[genome_build]} " \
			  "-protocol " + protocols + ' ' + \
			  "-operation " + operations + ' ' + \
			  "{config[annovar_db_storage]} "
			  "-remove ")
