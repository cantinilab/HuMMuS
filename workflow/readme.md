## HuMMuS workflow

In this directory, there is a snakemake workflow to run the full pipeline in R. You may modify the configuration adding seeds to perform the multilayer network analysis that you wish.
It requires two workflow parameters that you can see ann example in config/config.yaml:
1. **working_dir**: Path to the directory where the hummus files and directories will be stored
2. **configuration_file**: Path to the configuration file in json with the details of the execution and the experiment you wish to perform.

#### Configuration
The configuration file contains the following information:
1. **resources**: a dictionary containing the number of cores you want to assign when performing the analytical functions. Example: `{ "cores": 4 }`
2. **assays_rda_path**: a path to the Seurat assays .rda file
3. **network_methods**: The method that you want to use to obtain each multiplex network for TF, gene and peaks. The default values are in the config example. Example: <code>
    {
        "tf": "Omnipath", 
        "gene": "GENIE3", 
        "peak": "cicero"
    }
</code>

4. **analysis**: You may include here the inferred types of networks that you want to retrieve from the multilayer network analysis. You may also configure subset lists of genes and/or transcription factors to filter the results. You can also choose a custom name for the report. The output files of the analysis will be saved in working_dir/reports.
    Example: <code>
    {
        "grn": {
            "seeds": {
                "gene": "all", 
                "tf": "all"
            },
            "filename": "ranked_grn.tsv"
        }, 
        "target_gene": {
            "seeds": {
                "gene": "all", 
                "tf": "all"
            },
            "filename": "ranked_gene.tsv"
        }, 
        "enhancer": {
            "seeds": {
                "gene": "all",
                "tf": "all"
            },
            "filename": "ranked_enhancer.tsv"
        },
        "binding_region": {
            "seeds": {
                "gene": "all", 
                "tf": "all"
            },
            "filename": "ranked_binding_region.tsv"
        }
    }
    </code>

#### Usage:
- In the folder container/ you will find instructions to build a singularity image with all needed dependencies for HuMMuS work correctly.
- Run command example attaching a directory in the singularity image: `snakemake --cores 4 --use-singularity --singularity-args "--bind  '/mnt/data_volume/home/user//:/mnt/data_volume/home/user/'"`
