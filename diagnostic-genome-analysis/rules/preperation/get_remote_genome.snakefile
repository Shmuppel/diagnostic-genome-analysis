from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_remote_genome:
    """
    Retrieves a reference genome through HTTP transfer.

    Input: 
        HTTP RemoteProvider object with remote path to reference genome,
        taken from config.yaml.

    Output: 
        Reference genome in Fasta format.
    """
    input:
        HTTP.remote(config["remote_genome_path"], keep_local=True)
    output:
        config["genome"]
    group: "preperation"
    run:
        output_path = config["genome"]
        shell("mv {input} {output_path}")