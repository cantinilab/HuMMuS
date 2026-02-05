import hummuspy as hummus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

cpath = "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/hummus_cantinilab_issues/data/config.yml"
config = hummus.config.open_config(cpath)

multilayer_f = "./test_multiplex_genes/" # the general folder all the files
config_name = "extended_config.yaml" # the beautiful name you chose for your config
list_seeds = [] # Seeds from which we start the exploration. Generally all the nodes from the startin layer, indicated in 'eta' (e.g.: TFs) 

df = hummus.core_grn.define_grn_from_config(
        multilayer_f,
        config,
        "all",
        save=False, # Do we want to save the results on disk
        output_f=None, # Name of the result file IF save on disk
        return_df=True, # return in console the results
                                       # (e.g.: If you're interested only in TF-genes interactions and not in peaks-related scores, put ['RNA']
        njobs=1 # How many cores do you wanna use ?
)

# adding network

config['multiplex']['drugs'] = {
    'layers':  ["multiplex/drug_drug_network.tsv"],
    'graph_type': "01"
}
config['bipartite']["rna_drugs.tsv"] = {
    'source': "RNA",
    'target': "drugs",
    'graph_type': "00"
}

# Exploring network

genes = pd.read_csv("./data/multiplex/RNA/RNA_GENIE3.tsv", sep="\t", header=None)
genes = list(genes[1].unique())
config = hummus.config.process_config(config, multilayer_folder = multilayer_f )
config['seeds']+=genes # list(TFs) #["ARID3A", "ARID3B", "ARID2", "CLOCK", "ARNT2"]

df = hummus.explore_network.compute_multiple_RandomWalk( **config, save=False, n_jobs=4 )
