import typing
import os
import yaml
import re
import numpy as np
import pandas as pd
import multixrank as mxr
import time
from joblib import Parallel, delayed
from tqdm import tqdm

def make_values_list(values, types=(str, bool, int, float)):
    """Transform layer type to be sure a list of values is returned."""
    if type(values) == list:
        return values
    elif type(values)==dict:
        return list(values.values())
    elif type(values) in types:
        return [values]
    else:
        raise TypeError('Layer name(s) {} should be given through a {},'.format(values, types)+\
                        'a list, or a dictionary(where the key would be the file_paths)')
    
def define_minimal_config(
    tf_layers: typing.Union[str, list[str], dict[str]],
    atac_layers: typing.Union[str, list[str], dict[str]],
    rna_layers: typing.Union[str, list[str], dict[str]],
    tf_atac_links: str, 
    atac_rna_links: str,
    seed_path: str = 'seeds/seeds.txt',
    folder_tf_layers = 'layer_TFS',
    folder_atac_layers = 'layer_PEAKS',
    folder_rna_layers = 'layer_GENES',
    tf_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True, 
    rna_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True,
    tf_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    rna_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False,
    tf_atac_links_weighted: bool = False,
    atac_rna_links_weighted: bool = False,
    tf_atac_links_directed: bool = False, 
    atac_rna_links_directed: bool = False,
    r:float = 0.7,
    self_loops = 0):
    
    """eta and lambda will be added by specific tasks"""
    tf_layers = make_values_list(tf_layers)
    atac_layers = make_values_list(atac_layers)
    rna_layers = make_values_list(rna_layers)
    
    tf_layers_weighted = make_values_list(tf_layers_weighted)
    atac_layers_weighted = make_values_list(atac_layers_weighted)
    rna_layers_weighted = make_values_list(rna_layers_weighted)
    
    tf_layers_directed = make_values_list(tf_layers_directed)
    atac_layers_directed = make_values_list(atac_layers_directed)
    rna_layers_directed = make_values_list(rna_layers_directed)
    
    len_assert_msg = "number of layers not equal to number of graph is_weighted arguments: "
    assert len(tf_layers) == len(tf_layers_weighted), len_assert_msg+'{}, {}'.format(tf_layers, tf_layers_weighted)
    assert len(atac_layers) == len(atac_layers_weighted),  len_assert_msg+'{}, {}'.format(atac_layers, atac_layers_weighted)
    assert len(rna_layers) == len(rna_layers_weighted),  len_assert_msg+'{}, {}'.format(rna_layers, rna_layers_weighted)

    len_assert_msg = "number of layers not equal to number of graph is_directed arguments: "
    assert len(tf_layers) == len(tf_layers_directed), len_assert_msg+'{}, {}'.format(tf_layers, tf_layers_directed)
    assert len(atac_layers) == len(atac_layers_directed),  len_assert_msg+'{}, {}'.format(atac_layers, atac_layers_directed)
    assert len(rna_layers) == len(rna_layers_directed),  len_assert_msg+'{}, {}'.format(rna_layers, rna_layers_directed)                  

    
    assert tf_atac_links != atac_rna_links,  len_assert_msg+'The two bipartites {}, {} cannot have the same name'.format(
        tf_atac_links, 
        atac_rna_links)                  
                                                                                    

    tf_layers_types = [str(int(tf_layers_directed[i]))\
                       +str(int(tf_layers_weighted[i]))\
                              for i in range(len(tf_layers))]

    atac_layers_types = [str(int(atac_layers_directed[i]))\
                         +str(int(atac_layers_weighted[i]))\
                              for i in range(len(atac_layers))]
                                                                                      
    rna_layers_types = [str(int(rna_layers_directed[i]))\
                        +str(int(rna_layers_weighted[i]))\
                              for i in range(len(rna_layers))]
    config = dict()
    config['multiplex'] = dict()
    config['multiplex'][folder_tf_layers] = {
        'layers': [folder_tf_layers+'/'+tf_layer for tf_layer in tf_layers],
        'graph_type': tf_layers_types}

    config['multiplex'][folder_atac_layers] = {
        'layers': [folder_atac_layers+'/'+atac_layer for atac_layer in atac_layers],
        'graph_type': atac_layers_types}

    config['multiplex'][folder_rna_layers] = {
        'layers': [folder_rna_layers+'/'+rna_layer for rna_layer in rna_layers],
        'graph_type': rna_layers_types}
    # add delta and tau parameters

    tf_atac_links_type = str(int(tf_atac_links_directed))\
        + str(int(tf_atac_links_weighted))

    atac_rna_links_type = str(int(atac_rna_links_directed))\
        + str(int(atac_rna_links_weighted))

    config['bipartite'] = {}
    config['bipartite'][tf_atac_links] = {'source': folder_tf_layers,
                                          'target': folder_atac_layers,
                                          'graph_type': tf_atac_links_type}

    config['bipartite'][atac_rna_links] = {'source': folder_atac_layers,
                                           'target': folder_rna_layers,
                                           'graph_type': atac_rna_links_type}

    config['seed'] = seed_path
    config['r'] = r
    config['self_loops'] = self_loops

    return config



def detailed_config(
    request: str,
    filename: str,
    tf_layers: typing.Union[str, list[str], dict[str]],
    atac_layers: typing.Union[str, list[str], dict[str]],
    rna_layers: typing.Union[str, list[str], dict[str]],
    tf_atac_links: str, 
    atac_rna_links: str,
    seed_path: str = 'seeds/seeds.txt',
    folder_tf_layers = 'layer_TFS',
    folder_atac_layers = 'layer_PEAKS',
    folder_rna_layers = 'layer_GENES',
    tf_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True, 
    rna_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True,
    tf_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    rna_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False,
    tf_atac_links_weighted: bool = False,
    atac_rna_links_weighted: bool = False,
    tf_atac_links_directed: bool = False, 
    atac_rna_links_directed: bool = False,
    r:float = 0.7,
    self_loops = 0, 
    return_config = False):
    
    config = define_minimal_config(
        tf_layers=tf_layers,
        atac_layers=atac_layers,
        rna_layers=rna_layers,
        tf_atac_links=tf_atac_links, 
        atac_rna_links=atac_rna_links,
        seed_path=seed_path,
        folder_tf_layers = folder_tf_layers,
        folder_atac_layers = folder_atac_layers,
        folder_rna_layers = folder_rna_layers,
        tf_layers_weighted=tf_layers_weighted, 
        atac_layers_weighted=atac_layers_weighted, 
        rna_layers_weighted=rna_layers_weighted,
        tf_layers_directed=tf_layers_directed, 
        atac_layers_directed=atac_layers_directed, 
        rna_layers_directed=rna_layers_directed,
        tf_atac_links_weighted=tf_atac_links_weighted,
        atac_rna_links_weighted=atac_rna_links_weighted,
        tf_atac_links_directed=tf_atac_links_directed, 
        atac_rna_links_directed=atac_rna_links_directed,
        r=r,
        self_loops = self_loops)
    
    if request.upper()=='GRN':
        layers_eta = [1, 0, 0]
        
        #Check direction of lamb entry in the last multixrank version !!!
        go_to_tfs   = [0, '1/3', 0]
        go_to_peaks = [1, '1/3', 0.5]
        go_to_genes = [0, '1/3', 0.5]

    elif request.enhancers()=='enhancers':
        layers_eta = [0, 0, 1]

        #Check direction of lamb entry in the last multixrank version !!!
        go_to_tfs   = [0, 0, 0]
        go_to_peaks = [1, 1, 1]
        go_to_genes = [0, 0, 0]

    elif request.enhancers()=='target_regions':
        layers_eta = [1, 0, 0]

        #Check direction of lamb entry in the last multixrank version !!!
        go_to_tfs   = [0, 0, 0]
        go_to_peaks = [1, 1, 1]
        go_to_genes = [0, 0, 0]

    layers_proba = [go_to_tfs, go_to_peaks, go_to_genes]
    layers_names = [folder_tf_layers, folder_atac_layers, folder_rna_layers]
    
    #sort transition matrix based on order of layers in the config file, aka alphabetical order
    config['lamb'] = [[val for _, val in sorted(zip(layers_names, go_to_proba))]\
                      for _, go_to_proba in sorted(zip(layers_names, layers_proba))]
    config['eta']  = [eta for _, eta  in sorted(zip(layers_names, layers_eta))]
    
    with open(filename, 'w') as f:
        yaml.dump(config, f)

    if return_config:
        return config


def compute_multiple_RandomWalk(
    multilayer_f,
    config_name,
    output_f,
    list_seeds,
    config_folder='config',
    save=True,
    Return=False,
    spec_layer_result_saved='all',
    njobs=1):

    ranking_all_dfs = pd.DataFrame(columns = ['layer', 'target', 'path_layer', 'score', 'tf'])

    l_ranking_df = Parallel(n_jobs=njobs)(delayed(compute_RandomWalk)(multilayer_f, config_name, seeds, config_folder, spec_layer_result_saved) for seeds in tqdm(list_seeds))
    ranking_all_dfs = pd.concat([ranking_all_dfs]+l_ranking_df)

    if save:
        ranking_all_dfs.sort_values(by='score').to_csv(output_f, sep='\t', index=False, header=True)
    if Return:
        return ranking_all_dfs

def compute_RandomWalk(
    multilayer_f,
    config_name,
    seeds,
    config_folder='config',
    spec_layer_result_saved='all',
    unnamed=False,
    njobs=1):

    # seeds file names
    seeds = make_values_list(seeds)
    if unnamed is True:
        seeds_filename = 'seeds.txt'
        if njobs>1:
            raise Exception("Impossible to use unnamed seeds files while parallelising random walks.")
    else:
        seeds_filename = '_'.join(seeds)

    # write seeds file
    with open(multilayer_f+'/seeds/'+seeds_filename+'.txt', 'w') as f:
        f.write('\n'.join(seeds)+'\n')        
    
    # config file personalised with seed file
    with open(multilayer_f+'/{}/'.format(config_folder)+config_name, 'r') as f:
        config = yaml.load(f, Loader=yaml.BaseLoader)
        config['seed'] = 'seeds/'+seeds_filename+'.txt'
    with open(multilayer_f+'/{}/'.format(config_folder)+seeds_filename+'_'+config_name, 'w') as f:
        yaml.dump(config, f)
        
    # multixrank
    multixrank_obj = mxr.Multixrank(config=multilayer_f+'/{}/'.format(config_folder)+seeds_filename+'_'+config_name, wdir=multilayer_f)
    ranking_df = multixrank_obj.random_walk_rank()
    
    # and filter df results andadd seeds name
    ranking_df['tf'] = '_'.join(seeds)
    ranking_df = ranking_df[ranking_df.score > 0]  # ??
    ranking_df.columns = ['layer', 'target', 'path_layer', 'score', 'seed']
    if spec_layer_result_saved != 'all':
        if type(spec_layer_result_saved)==str:
            spec_layer_result_saved = [spec_layer_result_saved]
        ranking_df = ranking_df[ranking_df['layer'].isin(spec_layer_result_saved)]

    return ranking_df


def define_grn(
    multilayer_f: str,
    TFs: typing.Union[str, list[str], dict[str]] = 'all',
    save: bool = True,
    output_f_grn = 'RW_grn.tsv',
    Return: bool = False,
    njobs:int = 1,
    target_saved: typing.Union['genes', 'peaks', 'tfs'] = 'genes',
    tf_layers: typing.Union[str, list[str], dict[str]] = 'tf_edges.tsv',
    atac_layers: typing.Union[str, list[str], dict[str]] = 'peak_edges.tsv',
    rna_layers: typing.Union[str, list[str], dict[str]] = 'gene_edges.tsv',
    tf_atac_links: str = 'bipartite/tfs2peaks.tsv', 
    atac_rna_links: str = 'bipartite/peaks2genes.tsv',
    config_name='grn_config.yml',
    config_folder='config',
    seed_path: str = 'seeds/seeds.txt',
    folder_tf_layers: str = 'multiplex/layer_TFS',
    folder_atac_layers: str = 'multiplex/layer_PEAKS',
    folder_rna_layers: str = 'multiplex/layer_GENES',
    tf_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True, 
    rna_layers_weighted: typing.Union[bool, list[bool], dict[bool]] = True,
    tf_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    atac_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False, 
    rna_layers_directed: typing.Union[bool, list[bool], dict[bool]] = False,
    tf_atac_links_weighted: bool = False,
    atac_rna_links_weighted: bool = False,
    tf_atac_links_directed: bool = False, 
    atac_rna_links_directed: bool = False,
    r:float = 0.7,
    self_loops: bool = 0):

    detailed_config(
        request='grn',
        filename=multilayer_f+'/'+config_folder+'/'+config_name,
        tf_layers=tf_layers,
        atac_layers=atac_layers,
        rna_layers=rna_layers,
        tf_atac_links=tf_atac_links,
        atac_rna_links=atac_rna_links,
        seed_path=seed_path,
        folder_tf_layers=folder_tf_layers,
        folder_atac_layers=folder_atac_layers,
        folder_rna_layers=folder_rna_layers,
        tf_layers_weighted=tf_layers_weighted,
        atac_layers_weighted=atac_layers_weighted,
        rna_layers_weighted=rna_layers_weighted,
        tf_layers_directed=tf_layers_directed,
        atac_layers_directed=atac_layers_directed,
        rna_layers_directed=rna_layers_directed,
        tf_atac_links_weighted=tf_atac_links_weighted,
        atac_rna_links_weighted=atac_rna_links_weighted,
        tf_atac_links_directed=tf_atac_links_directed,
        atac_rna_links_directed=atac_rna_links_directed,
        r=r,
        self_loops=self_loops,
        return_config=False)
    
    if TFs == 'all':
        TFs = list(pd.read_csv(multilayer_f + '/' + tf_atac_links, sep="\t", header=None)[0].str.strip().unique())
    TFs = make_values_list(TFs)
    
    rna_layers = make_values_list(rna_layers)
    print([folder_rna_layers+'/'+rna_layer for rna_layer in rna_layers])

    if Return is True:
        grn = compute_multiple_RandomWalk(
            multilayer_f=multilayer_f,
            config_name=config_name,
            output_f=output_f_grn,
            list_seeds=TFs,
            config_folder=config_folder,
            save=save,
            Return=Return,
            spec_layer_result_saved=folder_rna_layers,
            njobs=njobs)
        return grn
    else:
        grn = compute_multiple_RandomWalk(
            multilayer_f=multilayer_f,
            config_name=config_name,
            output_f=output_f_grn,
            list_seeds=TFs,
            config_folder=config_folder,
            save=save,
            Return=Return,
            spec_layer_result_saved=folder_rna_layers,
            njobs=njobs)


# multilayer_f = '../../flattened_networks/ML_hESC_Chen_GeneNW_all_Peaks2Genes_UP0.5K_DOWN0.5K_nofilt_nofilt_TFlayer_nolinks'
# define_grn(multilayer_f, njobs=45)
