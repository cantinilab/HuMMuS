# import typing
# import os
# import time
import yaml
import numpy as np
import pandas as pd
import multixrank as mxr
from joblib import Parallel, delayed
from tqdm import tqdm
# from typing import Union
import hummuspy.config


def compute_multiple_RandomWalk(
        multilayer_f,
        config_name,
        list_seeds,
        config_folder='config',
        save=True,
        output_f=None,
        return_df=False,
        spec_layer_result_saved='all',
        njobs=1):
    """Compute random walks for a list of seeds.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config_name : str
        Name of the config file.
    output_f : str
        Name of the output file.
    list_seeds : list
        List of seeds.
    config_folder : str, optional
        Name of the config folder. The default is 'config'.
    save : bool, optional
        Save the result. The default is True.
    return_df : bool, optional
        Return the result. The default is False.
    spec_layer_result_saved : str, optional
        Name of the layer to save. The default is 'all'.
        Specify here if you want to keep only the scores of going to one of the layer 
            (e.g.: If you're interested only in TF-genes 
            interactions and not in peaks-related scores, put ['RNA']
    njobs : int, optional
        Number of jobs. The default is 1.

    Returns
    -------
    ranking_all_dfs : pd.DataFrame
        Dataframe containing the result of the random walks.
        Structure:
            layer : str
                Name of the target layer.
            target : str
                Name of the target.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            seed : str
                Name of the seed.

    Examples
    --------
    >>> import hummuspy
    >>> multilayer_f = 'path/to/multilayer/folder'
    >>> config_folder = 'config'
    >>> config_name = 'hummuspy.config.yml'
    >>> list_seeds = ['seed1', 'seed2']
    >>> df = compute_multiple_RandomWalk(multilayer_f,
                                    config_name,
                                    list_seeds,
                                    config_folder=config_folder,
                                    save=False,
                                    return_df=True,
                                    spec_layer_result_saved='all', # or 'TF'
                                    njobs=5)
    """

    ranking_all_dfs = pd.DataFrame(columns=['layer',
                                            'target',
                                            'path_layer',
                                            'score',
                                            'seed'])

    l_ranking_df = Parallel(n_jobs=njobs)(delayed(compute_RandomWalk)
                                          (multilayer_f=multilayer_f,
                                           config_name=config_name,
                                           seeds=seeds,
                                           config_folder=config_folder,
                                           save=False,
                                           spec_layer_result_saved=spec_layer_result_saved)
                                          for seeds in tqdm(list_seeds,
                                                            position=0,
                                                            leave=True))
    ranking_all_dfs = pd.concat([ranking_all_dfs]+l_ranking_df)
    ranking_all_dfs = ranking_all_dfs.sort_values(by='score', ascending=False)
    ranking_all_dfs = ranking_all_dfs.reset_index(drop=True)

    if save:
        assert output_f is not None, 'You need to provide an output_f name' +\
            ' to save the random walks result'
        ranking_all_dfs.to_csv(output_f, sep='\t', index=False, header=True)
    if return_df:
        return ranking_all_dfs


def compute_RandomWalk(
        multilayer_f,
        config_name,
        seeds,
        seeds_filename='auto',
        seeds_folder='seed',
        config_folder='config',
        save=True,
        output_f=None,
        return_df=True,
        spec_layer_result_saved='all',
        njobs=1):
    """Compute random walks for a list of seeds.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config_name : str
        Name of the config file.
    seeds : list
        List of seeds.
    config_folder : str, optional
        Name of the config folder. The default is 'config'.
    spec_layer_result_saved : str, optional
        Name of the layer to save. The default is 'all'.
    unnamed : bool, optional
        If True, the seeds file will be named 'seeds.txt'.
        The default is False.
    njobs : int, optional
        Number of jobs. The default is 1.

    Returns
    -------
    ranking_df : pd.DataFrame
        Dataframe containing the result of the random walk.
        Structure:
            layer : str
                Name of the target layer.
            target : str
                Name of the target.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            seed : str
                Name of the seed.

    Examples
    --------
    >>> import hummuspy
    >>> multilayer_f = 'path/to/multilayer/folder'
    >>> config_folder = 'config'
    >>> config_name = 'hummuspy.config.yml'
    >>> seed = 'seed1'
    >>> df = compute_RandomWalk(multilayer_f,
                                config_name,
                                seed,
                                # seeds_filename = 'auto'/'your_name.txt'
                                config_folder=config_folder,
                                spec_layer_result_saved='all', # or 'TF'
                                njobs=5)
    """
    seeds_filename
    # seeds file names
    seeds = hummuspy.config.make_values_list(seeds)
    if seeds_filename != 'auto':
        if njobs > 1:
            raise Exception("Impossible to use only one seeds filename while" +
                            " parallelising random walks." +
                            "\nTry seeds_filename = 'auto', or njobs=1.")
    else:
        seeds_filename = '_'.join(seeds)

    # write seeds file
    with open(multilayer_f + "/" + seeds_folder + "/" +
              seeds_filename + '.txt', 'w') as f:
        f.write('\n'.join(seeds)+'\n')

    # config file personalised with seed file
    with open(multilayer_f+'/{}/'.format(config_folder)+config_name, 'r') as f:
        config = yaml.load(f, Loader=yaml.BaseLoader)
        config['seed'] = seeds_folder+'/'+seeds_filename+'.txt'
    with open(multilayer_f+'/{}/'.format(config_folder)
              + seeds_filename + '_' + config_name, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    # multixrank
    multixrank_obj = mxr.Multixrank(
        config=multilayer_f + '/' + config_folder + '/'
        + seeds_filename + '_' + config_name,
        wdir=multilayer_f)
    ranking_df = multixrank_obj.random_walk_rank()

    # and filter df results and add seeds name
    ranking_df['seed'] = seeds_filename
    ranking_df = ranking_df[ranking_df.score > 0]  # ??
    ranking_df.columns = ['layer', 'target', 'path_layer', 'score', 'seed']
    if spec_layer_result_saved != 'all':
        if type(spec_layer_result_saved) == str:
            spec_layer_result_saved = [spec_layer_result_saved]
        ranking_df = ranking_df[ranking_df['layer'].isin(
            spec_layer_result_saved)]

    if save:
        assert output_f is not None, 'You need to provide an output_f name' +\
            ' to save the random walks result'
        ranking_df.to_csv(output_f, sep='\t', index=False, header=True)
    if return_df:
        return ranking_df


#############################################
# 1/4 Define GRN config and compute results #
#############################################
def define_grn_from_config(
        multilayer_f,
        config,
        gene_list=None,
        tf_list=None,
        config_name='grn_config.yml',
        config_folder='config',
        tf_multiplex: str = 'TF',
        peak_multiplex: str = 'peaks',
        rna_multiplex: str = 'RNA',
        update_config=True,
        save=False,
        return_df=True,
        output_f=None,
        njobs=1):
    """Define a GRN from a multilayer network and a config file.
    Random walks are computed for each gene in the gene list and we keep
    the probability to reach each TF in the TF list.
    You can provide a list of genes and TFs to restrict the GRN.
    The gene_list is used as individual seed for computing the random walks.
    The list of TFs is used after the random walks, filtering the results to
    only the TFs of interest.
    You can choose to save the result in a file and/or return it.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config : dict
        Config dictionnary.
    gene_list : list, optional
        List of genes. The default is 'all'.
    tf_list : list, optional
        List of TFs. The default is 'all'.
    config_name : str, optional
        Name of the config file that will be saved.
        The default is 'grn_config.yml'.
    config_folder : str, optional
        Name of the config folder where the config will be save.
        ! For each seed (sometimes thousands), a file should be created in this
        folder. The default is 'config'.
    tf_multiplex : str, optional
        Name of the TF multiplex. The default is 'TF'.
    peak_multiplex : str, optional
        Name of the peak multiplex. The default is 'peaks'.
    rna_multiplex : str, optional
        Name of the RNA multiplex. The default is 'RNA'.
    update_config : bool, optional
        Update the config file. The default is True ; if False, the config
        file won't be updated for the values of eta and lamb.
    save : bool, optional
        Save the result. The default is False. If True, you need to provide
        an output_f name to save the GRN result.
    return_df : bool, optional
        Return the result. The default is True.
    output_f : str, optional
        Name of the output file. The default is None. Only used if save=True.
    njobs : int, optional
        Number of jobs. The default is 1. If >1, the seeds will be saved in
        different files (in the multilayer subfolder 'seed') and the random
        walks will be parallelised.
    Returns
    -------
    df : pd.DataFrame
        Dataframe containing the random walks's results that defines the GRN.
        Columns:
            layer : str
                Name of the target layer.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            gene : str
                Name of the gene-seed.
            tf : str
                Name of the TF-target.

    """
    # store mutliplex already because it will be when saving yaml file,
    # while eta and lambda won't.
    config['multiplex'] = {k: config['multiplex'][k]
                           for k in sorted(config['multiplex'].keys())}

    if update_config:
        eta = hummuspy.config.get_single_layer_eta(config,
                                                   rna_multiplex)

        lamb = hummuspy.config.get_grn_lamb(config,
                                            tf_multiplex,
                                            peak_multiplex,
                                            rna_multiplex,
                                            draw=False)
        config['eta'] = eta
        config['lamb'] = lamb
        # config = hummuspy.config.setup_proba_config(config, eta, lamb)

    config_path = multilayer_f+'/'+config_folder+'/'+config_name
    hummuspy.config.save_config(config, config_path)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            gene_list = np.unique(np.concatenate([gene_list,
                                                  layer_nodes]))

    df = compute_multiple_RandomWalk(multilayer_f,
                                     config_name=config_name,
                                     output_f=output_f,
                                     list_seeds=gene_list,
                                     config_folder=config_folder,
                                     save=False,
                                     return_df=return_df,
                                     spec_layer_result_saved=tf_multiplex,
                                     njobs=njobs)

    df['gene'] = df['seed']
    df['tf'] = df['target']
    del df['target']
    del df['seed']

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            tf_list = np.unique(np.concatenate([tf_list,
                                                layer_nodes]))
        tf_list = tf_list[tf_list != 'fake_node']

    # Add normalisation ?
    df = df[df['tf'].isin(tf_list)]

    if save is True:
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the GRN result'
        df.sort_values(by='score', ascending=False).to_csv(output_f,
                                                           sep='\t',
                                                           index=False,
                                                           header=True)
    if return_df:
        return df


###############################################################################
#############################################
# 2/4 Define enhancers config and compute results #
#############################################
def define_enhancers_from_config(
        multilayer_f,
        config,
        gene_list=None,
        peak_list=None,
        config_name='enhancers_config.yml',
        config_folder='config',
        tf_multiplex: str = 'TF',
        peak_multiplex: str = 'peaks',
        rna_multiplex: str = 'RNA',
        update_config=True,
        save=False,
        return_df=True,
        output_f=None,
        njobs=1):
    """Return enhancers prediction from a multilayer network and a config file.
    Random walks are computed for each gene in the gene list and we keep
    the probability to reach each peak in the peak list.
    You can provide a peak_list and a gene_list to restrict the predictions.
    The gene_list is used as individual seed for computing the random walks.
    The list of peaks is used after the random walks, filtering the results to
    only the peaks of interest.
    You can choose to save the result in a file and/or return it.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config : dict
        Config dictionnary.
    gene_list : list, optional
        List of genes. The default is 'all'.
    peak_list : list, optional
        List of peaks. The default is 'all'.
    config_name : str, optional
        Name of the config file that will be saved.
        The default is 'enhancers_config.yml'.
    config_folder : str, optional
        Name of the config folder where the config will be save.
        ! For each seed (sometimes thousands), a file should be created in this
        folder. The default is 'config'.
    tf_multiplex : str, optional
        Name of the TF multiplex. The default is 'TF'.
    peak_multiplex : str, optional
        Name of the peak multiplex. The default is 'peaks'.
    rna_multiplex : str, optional
        Name of the RNA multiplex. The default is 'RNA'.
    update_config : bool, optional
        Update the config file. The default is True ; if False, the config
        file won't be updated for the values of eta and lamb.
    save : bool, optional
        Save the result. The default is False. If True, you need to provide
        an output_f name to save the predictions.
    return_df : bool, optional
        Return the result. The default is True.
    output_f : str, optional
        Name of the output file. The default is None. Only used if save=True.
    njobs : int, optional
        Number of jobs. The default is 1. If >1, the seeds will be saved in
        different files (in the multilayer subfolder 'seed') and the random
        walks will be parallelised.
    Returns
    -------
    df : pd.DataFrame
        Dataframe of the random walks's results that defines the predictions.
        Columns:
            layer : str
                Name of the target layer.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            gene : str
                Name of the gene-seed.
            peak : str
                Name of the peak-target.
    """
    # store mutliplex already because it will be when saving yaml file,
    # while eta and lambda won't.
    config['multiplex'] = {k: config['multiplex'][k]
                           for k in sorted(config['multiplex'].keys())}

    if update_config:
        # Indicate layer where to start the random walks : rna_multiplex
        eta = hummuspy.config.get_single_layer_eta(config,
                                                   rna_multiplex)

        # Define proba matrix to jump between layer : rna <--> peaks
        lamb = hummuspy.config.get_enhancers_lamb(config,
                                                  tf_multiplex,
                                                  peak_multiplex,
                                                  rna_multiplex,
                                                  draw=False)

        config['eta'] = eta
        config['lamb'] = lamb
        # config = hummuspy.config.setup_proba_config(config, eta, lamb)

    config_path = multilayer_f+'/'+config_folder+'/'+config_name
    hummuspy.config.save_config(config, config_path)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            gene_list = np.unique(np.concatenate([gene_list,
                                                  layer_nodes]))

    df = compute_multiple_RandomWalk(multilayer_f,
                                     config_name=config_name,
                                     output_f=output_f,
                                     list_seeds=gene_list,
                                     config_folder=config_folder,
                                     save=False,
                                     return_df=return_df,
                                     spec_layer_result_saved=peak_multiplex,
                                     # save only peaks proba
                                     njobs=njobs)

    df['gene'] = df['seed']
    df['peak'] = df['target']
    del df['target']
    del df['seed']

    if peak_list is None:
        peak_list = []
        for layer in config['multiplex'][peak_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            peak_list = np.unique(np.concatenate([peak_list,
                                                  layer_nodes]))

    # Add normalisation ?
    df = df[df['peak'].isin(peak_list)]

    if save is True:
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the enhancers prediction result.'
        df.sort_values(by='score', ascending=False).to_csv(output_f,
                                                           sep='\t',
                                                           index=False,
                                                           header=True)
    if return_df:
        return df


#########################################################
# 3/4 Define binding regions config and compute results #
#########################################################
def define_binding_regions_from_config(
        multilayer_f,
        config,
        tf_list=None,
        peak_list=None,
        config_name='binding_regions_config.yml',
        config_folder='config',
        tf_multiplex: str = 'TF',
        peak_multiplex: str = 'peaks',
        rna_multiplex: str = 'RNA',
        update_config=True,
        save=False,
        return_df=True,
        output_f=None,
        njobs=1):
    """Return binding regions prediction from a multilayer network and a config
    file. Random walks are computed for each TF in the TF list and we keep the
    probability to reach each peak in the peak list.
    You can provide a list of peaks and a tf_list to restrict the predictions.
    The list of TFs is used as individual seed for computing the random walks.
    The list of peaks is used after the random walks, filtering the results to
    only the peaks of interest.
    You can choose to save the result in a file and/or return it.


    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config : dict
        Config dictionnary.
    tf_list : list, optional
        List of TFs. The default is 'all'.
    peak_list : list, optional
        List of peaks. The default is 'all'.
    config_name : str, optional
        Name of the config file that will be saved.
        The default is 'binding_regions_config.yml'.
    config_folder : str, optional
        Name of the config folder where the config will be save.
        ! For each seed (sometimes thousands), a file should be created in this
         folder. The default is 'config'.
    tf_multiplex : str, optional
        Name of the TF multiplex. The default is 'TF'.
    peak_multiplex : str, optional
       Name of the peak multiplex. The default is 'peaks'.
    rna_multiplex : str, optional
        Name of the RNA multiplex. The default is 'RNA'.
    update_config : bool, optional
        Update the config file. The default is True ; if False, the config
        file won't be updated for the values of eta and lamb.
    save : bool, optional
        Save the result. The default is False. If True, you need to provide
        an output_f name to save the predictions.
    return_df : bool, optional
        Return the result. The default is True.
    output_f : str, optional
        Name of the output file. The default is None. Only used if save=True.
    njobs : int, optional
        Number of jobs. The default is 1. If >1, the seeds will be saved in
        different files (in the multilayer subfolder 'seed') and the random
        walks will be parallelised.

    Returns
    -------
    df : pd.DataFrame
        Dataframe of the random walks's results that defines the predictions.
        Columns:
            layer : str
                Name of the target layer.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            tf : str
                Name of the TF-seed.
            peak : str
                Name of the peak-target.
    """
    # store mutliplex already because it will be when saving yaml file,
    # while eta and lambda won't.
    config['multiplex'] = {k: config['multiplex'][k]
                           for k in sorted(config['multiplex'].keys())}

    if update_config:
        # Indicate layer where to start the random walks : rna_multiplex
        eta = hummuspy.config.get_single_layer_eta(config,
                                                   tf_multiplex)

        # Define proba matrix to jump between layer : rna <--> peaks
        lamb = hummuspy.config.get_binding_regions_lamb(config,
                                                        tf_multiplex,
                                                        peak_multiplex,
                                                        rna_multiplex,
                                                        draw=False)
        config['eta'] = eta
        config['lamb'] = lamb
        # config = hummuspy.config.setup_proba_config(config, eta, lamb)

    config_path = multilayer_f+'/'+config_folder+'/'+config_name
    hummuspy.config.save_config(config, config_path)

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            tf_list = np.unique(np.concatenate([tf_list,
                                                layer_nodes]))
        tf_list = tf_list[tf_list != 'fake_node']

    df = compute_multiple_RandomWalk(multilayer_f,
                                     config_name=config_name,
                                     output_f=output_f,
                                     list_seeds=tf_list,
                                     config_folder=config_folder,
                                     save=False,
                                     return_df=return_df,
                                     spec_layer_result_saved=peak_multiplex,
                                     # save only peaks proba
                                     njobs=njobs)

    df['tf'] = df['seed']
    df['peak'] = df['target']
    del df['target']
    del df['seed']

    if peak_list is None:
        peak_list = []
        for layer in config['multiplex'][peak_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            peak_list = np.unique(np.concatenate([peak_list,
                                                  layer_nodes]))

    # Add normalisation ?
    df = df[df['peak'].isin(peak_list)]

    if save is True:
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the enhancers prediction result.'
        df.sort_values(by='score', ascending=False).to_csv(output_f,
                                                           sep='\t',
                                                           index=False,
                                                           header=True)
    if return_df:
        return df


######################################################
# 4/4 Define target genes config and compute results #
######################################################
def define_target_genes_from_config(
        multilayer_f,
        config,
        gene_list=None,
        tf_list=None,
        config_name='target_genes_config.yml',
        config_folder='config',
        tf_multiplex: str = 'TF',
        peak_multiplex: str = 'peaks',
        rna_multiplex: str = 'RNA',
        update_config=True,
        save=False,
        return_df=True,
        output_f=None,
        njobs=1):
    """Return target genes prediction from a multilayer network and a config
    file. Random walks are computed for each TF in the TF list and we keep the
    probability to reach each gene in the gene list.
    You can provide a list of genes and a tf_list to restrict the predictions.
    The list of TFs is used as individual seed for computing the random walks.
    The list of genes is used after the random walks, filtering the results to
    only the genes of interest.
    You can choose to save the result in a file and/or return it.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config : dict
        Config dictionnary.
    gene_list : list, optional
        List of genes. The default is 'all'.
    tf_list : list, optional
        List of TFs. The default is 'all'.
    config_name : str, optional
        Name of the config file that will be saved.
        The default is 'target_genes_config.yml'.
    config_folder : str, optional
        Name of the config folder where the config will be save.
        ! For each seed (sometimes thousands), a file should be created in this
        folder. The default is 'config'.
    tf_multiplex : str, optional
        Name of the TF multiplex. The default is 'TF'.
    peak_multiplex : str, optional
        Name of the peak multiplex. The default is 'peaks'.
    rna_multiplex : str, optional
        Name of the RNA multiplex. The default is 'RNA'.
    update_config : bool, optional
        Update the config file. The default is True ; if False, the config
        file won't be updated for the values of eta and lamb.
    save : bool, optional
        Save the result. The default is False. If True, you need to provide
        an output_f name to save the predictions.
    return_df : bool, optional
        Return the result. The default is True.
    output_f : str, optional
        Name of the output file. The default is None. Only used if save=True.
    njobs : int, optional
        Number of jobs. The default is 1. If >1, the seeds will be saved in
        different files (in the multilayer subfolder 'seed') and the random
        walks will be parallelised.

    Returns
    -------
    df : pd.DataFrame
        Dataframe of the random walks's results that defines the predictions.
        Columns:
            layer : str
                Name of the target layer.
            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            tf : str
                Name of the TF-seed.
            gene : str
                Name of the gene-target.
    """
    # store mutliplex already because it will be when saving yaml file,
    # while eta and lambda won't.
    config['multiplex'] = {k: config['multiplex'][k]
                           for k in sorted(config['multiplex'].keys())}

    if update_config:
        eta = hummuspy.config.get_single_layer_eta(config,
                                                   tf_multiplex)

        lamb = hummuspy.config.get_target_genes_lamb(config,
                                                     tf_multiplex,
                                                     peak_multiplex,
                                                     rna_multiplex,
                                                     draw=False)

        config['eta'] = eta
        config['lamb'] = lamb
        # config = hummuspy.config.setup_proba_config(config, eta, lamb)

    config_path = multilayer_f+'/'+config_folder+'/'+config_name
    hummuspy.config.save_config(config, config_path)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            gene_list = np.unique(np.concatenate([gene_list,
                                                  layer_nodes]))

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pd.read_csv(multilayer_f+'/'+layer,
                                   sep='\t',
                                   header=None,
                                   index_col=None)

            layer_nodes = np.concatenate([np.unique(df_layer[0].values),
                                          np.unique(df_layer[1].values)])
            tf_list = np.unique(np.concatenate([tf_list,
                                                layer_nodes]))
        tf_list = tf_list[tf_list != 'fake_node']

    df = compute_multiple_RandomWalk(multilayer_f,
                                     config_name=config_name,
                                     output_f=output_f,
                                     list_seeds=tf_list,
                                     config_folder=config_folder,
                                     save=False,
                                     return_df=return_df,
                                     spec_layer_result_saved=rna_multiplex,
                                     njobs=njobs)

    df['tf'] = df['seed']
    df['gene'] = df['target']
    del df['target']
    del df['seed']

    # Add normalisation ?
    df = df[df['gene'].isin(gene_list)]

    if save is True:
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the GRN result'
        df.sort_values(by='score', ascending=False).to_csv(output_f,
                                                           sep='\t',
                                                           index=False,
                                                           header=True)
    if return_df:
        return df


def get_output_from_dicts(
        output_request: str,
        multilayer_f,
        multiplexes_list,
        bipartites_list,
        folder_multiplexes='multiplex',
        folder_bipartites='bipartite',
        gene_list=None,
        tf_list=None,
        peak_list=None,
        config_filename='config.yml',
        config_folder='config',
        tf_multiplex: str = 'TF',
        peak_multiplex: str = 'peaks',
        rna_multiplex: str = 'RNA',
        bipartites_type=('00', '00'),
        update_config=True,
        save=False,
        return_df=True,
        output_f=None,
        njobs=1):
    """
    Compute an output from a multilayer network and a config file, that can be
    chosen among ['grn', 'enhancers', 'binding_regions', 'target_genes'].

    It is a wrapper of the functions define_*_from_config, that are called
    depending on the output_request parameter.

    Parameters
    ----------
    output_request : ['grn', 'enhancers', 'binding_regions', 'target_genes']
        Type of output requested.
    multilayer_f : str
        Path to the multilayer folder.
    config : dict
        Config dictionnary.
    gene_list : list, optional
        List of genes. The default is 'all'.
    tf_list : list, optional
        List of TFs. The default is 'all'.
    config_name : str, optional
        Name of the config file. The default is 'config.yml'.
    config_folder : str, optional
        Name of the config folder. The default is 'config'.
    tf_multiplex : str, optional
        Name of the TF multiplex. The default is 'TF'.
    peak_multiplex : str, optional
        Name of the peak multiplex. The default is 'peaks'.
    rna_multiplex : str, optional
        Name of the RNA multiplex. The default is 'RNA'.
    update_config : bool, optional
        Update the config file. The default is True.
    save : bool, optional
        Save the result. The default is False.
    return_df : bool, optional
        Return the result. The default is True.
    output_f : str, optional
        Name of the output file. The default is None.
    njobs : int, optional
        Number of jobs. The default is 1.

    Returns
    -------ith open(self.config_path) as fin:
            self.config_dic = yaml.load(fin, Loader=yaml.BaseLoader)
    df : pd.DataFrame
        Dataframe containing the random walks's results that defines the GRN.
        Columns:
            layer : str
                Name of the target layer.

            path_layer : str
                Name of the layer of the path.
            score : float
                Score of the random walk.
            gene : str
                Name of the gene-seed.
            tf : str
                Name of the TF-target.

    """

    njobs = int(njobs)
    print('multiplexes_list : ', multiplexes_list)
    print('bipartites_list : ', bipartites_list)
    print('folder_multiplexes : ', folder_multiplexes)
    print('folder_bipartites : ', folder_bipartites)
    print('gene_list : ', gene_list)
    print('tf_list : ', tf_list)
    print('peak_list : ', peak_list)
    print('config_filename : ', config_filename)
    print('config_folder : ', config_folder)
    print('tf_multiplex : ', tf_multiplex)
    print('peak_multiplex : ', peak_multiplex)
    print('rna_multiplex : ', rna_multiplex)
    print('update_config : ', update_config)
    print('save : ', save)
    print('return_df : ', return_df)
    print('output_f : ', output_f)
    print('njobs : ', njobs)

    # Create general config file
    config = hummuspy.config.general_config(
        multiplexes=multiplexes_list,
        bipartites=bipartites_list,
        folder_multiplexes=folder_multiplexes,
        folder_bipartites=folder_bipartites,
        suffix='.tsv',
        self_loops=0,
        restart_prob=0.7,
        bipartites_type=bipartites_type,
        save_configfile=False,
        config_filename=config_filename)

    parameters = {
        'multilayer_f':   multilayer_f,
        'config':         config,
        'gene_list':      gene_list,
        'peak_list':      peak_list,
        'tf_list':        tf_list,
        'config_name':    config_filename,
        'config_folder':  config_folder,
        'tf_multiplex':   tf_multiplex,
        'peak_multiplex': peak_multiplex,
        'rna_multiplex':  rna_multiplex,
        'update_config':  update_config,
        'save':           save,
        'return_df':      return_df,
        'output_f':       output_f,
        'njobs':          njobs
    }

    if output_request == 'grn':
        del parameters['peak_list']
        df = define_grn_from_config(**parameters)

    elif output_request == 'enhancers':
        del parameters['tf_list']
        df = define_enhancers_from_config(**parameters)

    elif output_request == 'binding_regions':
        del parameters['gene_list']
        df = define_binding_regions_from_config(**parameters)

    elif output_request == 'target_genes':
        del parameters['peak_list']
        df = define_target_genes_from_config(**parameters)
    else:
        raise ValueError("Please select an output_request value in ('grn', 'enhancers', 'binding_regions', 'target_genes').")

    return df
