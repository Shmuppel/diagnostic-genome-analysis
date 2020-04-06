from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_remote_genome:
    """
    Retrieves a reference genome through HTTP transfer.

    Input: 
        HTTP RemoteProvider object with remote path to reference genome, taken from config.yaml.
    
    Output: 
        Reference genome in Fasta format.
        
    Threads:
        1, as this rule is a single download action it is single threaded. 
    """
    input:
        HTTP.remote(config["remote_genome_path"], keep_local=True)
    output:
        config["genome"]
    threads: 1
    priority: 1
    run:
        output_path = config["genome"]
        shell("mv {input} {output_path}")