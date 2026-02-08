import os
import hummuspy as hummus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

folder = "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/hummus_cantinilab_issues/HuMMuS/hummuspy/tests/data/"
tmp = os.path.join( folder, 'tmp_results')
dryrun_dir = os.path.join( folder, 'dryrun_files')

# ------------ test 1 [ok]
multiplex = {
    'RNA': {
        'graph_type': ['00'], 
        'layers': ['multiplex/RNA/GENIE3.tsv']
    },
    'peaks': {
        'graph_type': ['00'], 
        'layers': ['multiplex/TF/PPI.tsvmultiplex/peaks/Cicero.tsv']
    }
}
bipartite = {
    'bipartite/atac_rna.tsv': {
        'multiplex_right': 'RNA',
        'multiplex_left': 'peaks'
    }
}
bipartites_type = ['00']
config = hummuspy.config.general_config( multiplex, bipartite, bipartites_type = bipartites_type )

# ------------ test 2 [ok]
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

# ------------ test 3 [ok]
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

# adding network
config['multiplex']['drugs'] = {
    'layers':  ["./data/chen_multilayer/multiplex/drugs/drug_drug_network.tsv"],
    'graph_type': ["01"],
    'names': ['drugs_1']
}
config['bipartite']["bipartite/rna_drugs.tsv"] = {
    'edge_list_df': './data/chen_multilayer/bipartite/rna_drugs.tsv',
    'source': "RNA",
    'target': "drugs",
    'graph_type': "01"
}
config['eta'] = pd.Series([1, 0, 0, 0], index=['TF', 'peaks', 'RNA', 'drugs'])
config['lamb'] = hummus.config.initialise_lamb(config, value=0)

config['lamb'].loc["TF", "TF"] = 0.5 # probability to go to TF from TF
config['lamb'].loc["peaks", "TF"] = 0.5 # probability to go to peaks from TF

config['lamb'].loc["TF", "peaks"] = 1/3 # probability to go to TF FROM peaks
config['lamb'].loc["peaks", "peaks"] = 1/3 
config['lamb'].loc["RNA", "peaks"] = 1/3 

config['lamb'].loc["peaks", "RNA"] = 1/3  # probability to go to peaks FROM RNA
config['lamb'].loc["RNA", "RNA"] = 1/3 
config['lamb'].loc["drugs", "RNA"] = 1/3 

config['lamb'].loc["RNA", "drugs"] = 1/2
config['lamb'].loc["drugs", "drugs"] = 1/2

plpatha = os.path.join( tmp, "grn_plot.png" )
config['lamb'] = hummus.config.get_grn_lamb(config, draw=False, save_plot=True, name_plot=plpatha )
plpathb = os.path.join( tmp, "enhancer_plot.png" )
config['lamb'] = hummus.config.get_enhancers_lamb(config, draw=True, save_plot=True, name_plot=plpathb )
plpathc = os.path.join( tmp, "br_plot.png" )
config['lamb'] = hummus.config.get_binding_regions_lamb(config, draw=True, save_plot=True, name_plot=plpathc )
plpathd = os.path.join( tmp, "tg_plot.png" )
config['lamb'] = hummus.config.get_target_genes_lamb(config, draw=True, save_plot=True, name_plot=plpathd )

configExt_name = os.path.join( tmp, "extended_config.yaml") # the beautiful name you chose for your config
hummus.config.save_config(config, configExt_name )

assert os.path.isfile(configExt_name) == True
assert os.path.isfile(plpatha) == True
assert os.path.isfile(plpathb) == True
assert os.path.isfile(plpathc) == True
assert os.path.isfile(plpathd) == True

# ------------ test 4 [ok] 
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)
multilayer_folder = "./data/chen_multilayer"

save_flag = True
out_path = os.path.join( tmp, "ranked_grn_out.tsv" )
result = hummus.core_grn.define_grn_from_config(
        multilayer_folder,
        config,
        "all",
        save = save_flag, # Do we want to save the results on disk
        output_f = out_path, # Name of the result file IF save on disk
        return_df = True, # return in console the results
                                       # (e.g.: If you're interested only in TF-genes interactions and not in peaks-related scores, put ['RNA']
        return_config = True,
        njobs = 1 # How many cores do you wanna use ?
)
assert type(result) == tuple
df, config = result
assert len(df) == ?

# ------------ test 5
multilayer_folder = "./data/chen_multilayer"
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

save_flag = True
out_path = os.path.join( tmp, "ranked_enhancer_out.tsv" )
df, config = hummus.core_grn.define_enhancers_from_config(
        multilayer_folder,
        config,
        "all",
        save = save_flag, # Do we want to save the results on disk
        output_f = out_path, # Name of the result file IF save on disk
        return_df = True, # return in console the results
        return_config = True,
        njobs = 1 # How many cores do you wanna use ?
)
assert type(result) == tuple

# ------------ test 6
multilayer_folder = "./data/chen_multilayer"
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

save_flag = True
out_path = os.path.join( tmp, "ranked_bind_region_out.tsv" )
result = hummus.core_grn.define_binding_regions_from_config(
        multilayer_folder,
        config,
        "all",
        save = save_flag, # Do we want to save the results on disk
        output_f = out_path, # Name of the result file IF save on disk
        return_df = True, # return in console the results
                                       # (e.g.: If you're interested only in TF-genes interactions and not in peaks-related scores, put ['RNA']
        njobs = 1 # How many cores do you wanna use ?
)
assert type(result)==pandas.DataFrame

# ------------ test 7
multilayer_folder = "./data/chen_multilayer"
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

save_flag = True
out_path = os.path.join( tmp, "ranked_genes_out.tsv" )
df = hummus.core_grn.define_target_genes_from_config(
        multilayer_folder,
        config,
        "all",
        save = save_flag, # Do we want to save the results on disk
        output_f = out_path, # Name of the result file IF save on disk
        return_df = True, # return in console the results
                                       # (e.g.: If you're interested only in TF-genes interactions and not in peaks-related scores, put ['RNA']
        njobs = 1 # How many cores do you wanna use ?
)


# ------------ test 8
multilayer_folder = "./data/chen_multilayer"
cpath = os.path.join(folder, "chen_multilayer", "config", "grn_config.yml")
config = hummus.config.open_config(cpath)

# Exploring network
save_flag = True
out_path = os.path.join( tmp, "filtered_rw_out.tsv" )
rna_path = os.path.join( folder, 'chen_multilayer', 'multiplex', 'RNA', 'RNA_GENIE3.tsv' )
genes = pd.read_csv( rna_path, sep="\t", header=None)
genes = list(genes[1].unique())[:200]
config = hummus.config.process_config(config, multilayer_folder = multilayer_folder )
config['seeds'] = genes[:200]

df = hummus.explore_network.compute_multiple_RandomWalk( **config, save=save_flag, output_f=out_path, n_jobs=4 )
