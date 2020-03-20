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
    may be downloaded dynamically. Downloads are run in parallel using GNU parallel.
    
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
    run:
        command = "("
        for db in config["annovar_dbs"].keys():
            command += "{config[annovar_tool]}annotate_variation.pl " \
                       "-buildver {config[genome_build]} " \
                       "-downdb -webfrom annovar " + db + \
                       " {config[annovar_db_storage]} ; "

        command += ") | parallel"

        shell(command)
