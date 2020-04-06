def get_annovar_db_paths():
    """
    Retrieves the local paths where Annovar databases are/should be stored.

    Used so that databases only have to be downloaded once, and so that dependencies
    can be dynamically added / removed in the config.
    """
    # example db path: annovar databases folder + genome build prefix + annovar database file; for each database.
    return [config["annovar_db_storage"] + config["genome_build"] + '_' + config["annovar_dbs"][db]["file"]
            for db in config["annovar_dbs"]]

db_paths = get_annovar_db_paths()

rule get_remote_annovar_dbs:
    """
    Downloads Annovar variant annotation databases. 
    
    Databases that should be used are to be added in config.yaml, so that they
    may be downloaded dynamically. 
    
    Threads:
        1, with the assumption that the downloading of databases is usually going to be bandwidth or
        I/O limited rather than limited by computational power I've opted to use a Mutex/Semaphore. 
        This is enables the pipeline to reserve 1 thread to the download of databases whilst spending 
        the bulk of computational power on genome analysis tools. Feel free to increase the number of 
        threads if not bandwidth limited.
    
    Output:
        Database files as specified in config.yaml
    
    Run clarification: 
        # Run for every Annovar database in config.yaml
        <annovar tool path>annotate_variation.pl
        -buildver <genome build version>
        -downdb -webfrom annovar (download database) <database>
        <annovar databases storage folder>
    """
    output: db_paths
    log:
        "runs/pipeline_preperation/logs/get_remote_annovar_dbs.log"
    threads: min(config["reserve_annovar_db_thread"], workflow.cores)
    priority: 1
    run:
        # initiate a queue of databases to be downloaded.
        for db in config["annovar_dbs"].keys():
            command = "{config[annovar_tool]}annotate_variation.pl " \
                       "-buildver {config[genome_build]} " \
                       "-downdb -webfrom annovar " + db + " {config[annovar_db_storage]} "
            shell("(sem -j {threads} --id get_remote_annovar_dbs " + command + " ) >> {log}")

        # wait for all downloads to be completed before continuing.
        shell("(sem --wait --id get_remote_annovar_dbs ) >> {log}")
        # update local config, a thread no longer has to be reserved for this rule (will not change config.yaml).
        snakemake.utils.update_config(snakemake.workflow.config, {"reserve_annovar_db_thread": 0})
