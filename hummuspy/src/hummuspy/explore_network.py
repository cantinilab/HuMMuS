from hummuspy.create_multilayer import Multixrank
import os
import numpy
import pandas


def compute_RandomWalk(
        multiplex,
        bipartite,
        eta,
        lamb,
        seeds,
        self_loops=True,
        restart_proba=0.7,
        pr=None,
        save=True,
        output_f=None,
        return_df=True,
        spec_layer_result_saved='all',
        n_jobs=1):
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

    # multixrank
    multixrank_obj = Multixrank(
        multiplex=multiplex,
        bipartite=bipartite,
        eta=eta,
        lamb=lamb,
        seeds=seeds,
        self_loops=self_loops,
        restart_proba=restart_proba,
        pr=pr)
    ranking_df = multixrank_obj.random_walk_rank().sort_values(
        by='score',
        ascending=False)

    # and filter df results and add seeds name
    ranking_df['seed'] = '_'.join(seeds)
    ranking_df = ranking_df[ranking_df.score > 0]  # ??
    ranking_df.columns = ['layer', 'target', 'path_layer', 'score', 'seed']
    if spec_layer_result_saved != 'all':
        if type(spec_layer_result_saved) is str:
            spec_layer_result_saved = [spec_layer_result_saved]
        ranking_df = ranking_df[ranking_df['layer'].isin(
            spec_layer_result_saved)]

    if save:
        assert output_f is not None, 'You need to provide an output_f name' +\
            ' to save the random walks result'
        ranking_df.to_csv(output_f, sep='\t', index=False, header=True)
    if return_df:
        return ranking_df


def compute_multiple_RandomWalk(
        multiplex,
        bipartite,
        eta,
        lamb,
        seeds,
        self_loops=True,
        restart_proba=0.7,
        pr=None,
        save=True,
        output_f=None,
        return_df=True,
        spec_layer_result_saved='all',
        n_jobs=1):
    """Compute random walks for a list of seeds.

    Parameters
    ----------
    multilayer_f : str
        Path to the multilayer folder.
    config_name : strLINC01409
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

    # multixrank
    multixrank_obj = Multixrank(
        multiplex=multiplex,
        bipartite=bipartite,
        eta=eta,
        lamb=lamb,
        seeds=seeds,
        self_loops=self_loops,
        restart_proba=restart_proba,
        pr=pr)

    ranking_df = multixrank_obj.per_seed_random_walk_rank(
        n_jobs=n_jobs).sort_values(by='score', ascending=False)

    # and filter df results and add seeds name
    ranking_df = ranking_df[ranking_df.score > 0]  # ??
    ranking_df.columns = ['layer', 'target', 'path_layer', 'score', 'seed']
    if spec_layer_result_saved != 'all':
        if type(spec_layer_result_saved) is str:
            spec_layer_result_saved = [spec_layer_result_saved]
        ranking_df = ranking_df[ranking_df['layer'].isin(
            spec_layer_result_saved)]

    ranking_df = ranking_df.reset_index()

    if save:
        assert output_f is not None, 'You need to provide an output_f name' +\
            ' to save the random walks result'
        ranking_df.to_csv(output_f, sep='\t', index=False, header=True)
    if return_df:
        return ranking_df
