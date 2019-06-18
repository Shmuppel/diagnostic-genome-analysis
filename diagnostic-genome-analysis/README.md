# README #

This repository only contains the tools and rules (that is, set statements on how these tools are executed).
Certain samples and references are vastly too large to be stored in an online repository, and should be provided
manually. 

For this pipeline to work in its current state the following criteria should be met:

  * Samples that are to be subjected to analysis should be placed within resources/samples.
  * A reference genome should be provided under resources/reference.

Additionally a local Snakemake environment should be constructed around the Snakefile.
For more information on this see the [Snakemake documentation](http://snakemake.readthedocs.io/en/stable/).

#### Currently working on:

  * Incorporating more tools from the course guide into the pipeline. These tools include Varscan 2 
    and Annovar. Depending on the compatibility of subsequent tools some tools might be omitted to save time 
    and clarity
  * Creating a config file in either JSON or YAML to easily change variables (such as which variables to 
    analyse, or which reference genome to use) without editing the Snakefile.
  * Documentation of the Snakefile, its rules and the shell commands used to execute tools.