import numpy
from typing import Union
import pandas
import hummuspy.config
import hummuspy.explore_network
import hummuspy


def format_multilayer(
    TF_layer: Union[list, pandas.DataFrame],
    ATAC_layer: Union[list, pandas.DataFrame],
    RNA_layer: Union[list, pandas.DataFrame],
    TF_ATAC_bipartite: pandas.DataFrame,
    ATAC_RNA_bipartite: pandas.DataFrame,
    TF_layer_graph_type: Union[str, list] = '00',
    ATAC_layer_graph_type: Union[str, list] = '01',
    RNA_layer_graph_type: Union[str, list] = '01',
    TF_ATAC_bipartite_graph_type: str = '00',
    ATAC_RNA_bipartite_graph_type: str = '00'
):

    """ Format layers and bipartites data as needed for HuMMuS/MultiXrank

    ! The DataFrame should have 2-3 columns :
     ['source', 'target', 'weight'], or ['source', 'target'] !
    #TODO : explain graph_type !!
    Args:
        TF_layer (Union[list, pandas.DataFrame]): TF layer(s) edge list
        ATAC_layer (Union[list, pandas.DataFrame]): ATAC layer(s) edge list
        RNA_layer (Union[list, pandas.DataFrame]): RNA layer(s) edge list
        TF_ATAC_bipartite (pandas.DataFrame): TF-ATAC bipartite edge list
        ATAC_RNA_bipartite (pandas.DataFrame): ATAC-RNA bipartite edge list
        TF_layer_graph_type (Union[str, list], optional):
            Graph type for TF layer(s).
            Defaults to '00'.
        ATAC_layer_graph_type (Union[str, list], optional):
            Graph type for ATAC layer(s).
            Defaults to '01'.
        RNA_layer_graph_type (Union[str, list], optional):
            Graph type for RNA layer(s).
            Defaults to '01'.
        TF_ATAC_bipartite_graph_type (str, optional):
            Graph type for TF-ATAC bipartite.
            Defaults to '00'.
        ATAC_RNA_bipartite_graph_type (str, optional):
            Graph type for ATAC-RNA bipartite.
            Defaults to '00'.

    Returns:
        dict: Formatted data as input for MultiXrank-HuMMuS class
    """
    if type(TF_layer) is pandas.DataFrame:
        TF_layer = [TF_layer]
    if type(ATAC_layer) is pandas.DataFrame:
        ATAC_layer = [ATAC_layer]
    if type(RNA_layer) is pandas.DataFrame:
        RNA_layer = [RNA_layer]

    if type(TF_layer_graph_type) is str:
        TF_layer_graph_type = [TF_layer_graph_type]*len(TF_layer)
    if type(ATAC_layer_graph_type) is str:
        ATAC_layer_graph_type = [ATAC_layer_graph_type]*len(ATAC_layer)
    if type(RNA_layer_graph_type) is str:
        RNA_layer_graph_type = [RNA_layer_graph_type]*len(RNA_layer)

    assert len(TF_layer) == len(TF_layer_graph_type), \
        "TF_layer and TF_layer_graph_type must have the same length"
    assert len(ATAC_layer) == len(ATAC_layer_graph_type), \
        "ATAC_layer and ATAC_layer_graph_type must have the same length"
    assert len(RNA_layer) == len(RNA_layer_graph_type), \
        "RNA_layer and RNA_layer_graph_type must have the same length"

    keys_data = ['names', 'graph_type', 'layers']
    # Multiplex
    multiplex = {
        'TF': {key: [] for key in keys_data},
        'ATAC': {key: [] for key in keys_data},
        'RNA': {key: [] for key in keys_data}
    }

    for i, layer in enumerate(TF_layer):
        multiplex['TF']['names'].append(f'TF_{i}')
        multiplex['TF']['graph_type'].append(TF_layer_graph_type[i])
        multiplex['TF']['layers'].append(layer)

    for i, layer in enumerate(ATAC_layer):
        multiplex['ATAC']['names'].append(f'ATAC_{i}')
        multiplex['ATAC']['graph_type'].append(ATAC_layer_graph_type[i])
        multiplex['ATAC']['layers'].append(layer)

    for i, layer in enumerate(RNA_layer):
        multiplex['RNA']['names'].append(f'RNA_{i}')
        multiplex['RNA']['graph_type'].append(RNA_layer_graph_type[i])
        multiplex['RNA']['layers'].append(layer)

    # Bipartite
    bipartite = {
        'TF_ATAC': {
            'source': 'TF',
            'target': 'ATAC',
            'edge_list_df': TF_ATAC_bipartite,
            'graph_type': TF_ATAC_bipartite_graph_type},

        'ATAC_RNA': {
            'source': 'ATAC',
            'target': 'RNA',
            'edge_list_df': ATAC_RNA_bipartite,
            'graph_type': ATAC_RNA_bipartite_graph_type}
    }

    return multiplex, bipartite


def init_GRN(
    TF_layer: Union[list, pandas.DataFrame],
    ATAC_layer: Union[list, pandas.DataFrame],
    RNA_layer: Union[list, pandas.DataFrame],
    TF_ATAC_bipartite: pandas.DataFrame,
    ATAC_RNA_bipartite: pandas.DataFrame,
    TF_layer_graph_type: Union[str, list] = '00',
    ATAC_layer_graph_type: Union[str, list] = '10',
    RNA_layer_graph_type: Union[str, list] = '10',
    TF_ATAC_bipartite_graph_type: str = '00',
    ATAC_RNA_bipartite_graph_type: str = '00'
):

    multiplex, bipartite = format_multilayer(
        TF_layer=TF_layer,
        ATAC_layer=ATAC_layer,
        RNA_layer=RNA_layer,
        TF_ATAC_bipartite=TF_ATAC_bipartite,
        ATAC_RNA_bipartite=ATAC_RNA_bipartite,
        TF_layer_graph_type=TF_layer_graph_type,
        ATAC_layer_graph_type=ATAC_layer_graph_type,
        RNA_layer_graph_type=RNA_layer_graph_type,
        TF_ATAC_bipartite_graph_type=TF_ATAC_bipartite_graph_type,
        ATAC_RNA_bipartite_graph_type=ATAC_RNA_bipartite_graph_type
        )

    params = {
        'multiplex': multiplex,
        'bipartite': bipartite
        }

    lamb = hummuspy.config.get_grn_lamb(
        params,
        tf_multiplex='TF',
        peak_multiplex='ATAC',
        rna_multiplex='RNA')
    eta = hummuspy.config.get_single_layer_eta(
        params,
        'TF')

    return multiplex, bipartite, eta, lamb


def get_GRN(
    TF_layer: Union[list, pandas.DataFrame],
    ATAC_layer: Union[list, pandas.DataFrame],
    RNA_layer: Union[list, pandas.DataFrame],
    TF_ATAC_bipartite: pandas.DataFrame,
    ATAC_RNA_bipartite: pandas.DataFrame,
    seeds: Union[list, str],
    TF_layer_graph_type: Union[str, list] = '00',
    ATAC_layer_graph_type: Union[str, list] = '01',
    RNA_layer_graph_type: Union[str, list] = '01',
    TF_ATAC_bipartite_graph_type: str = '00',
    ATAC_RNA_bipartite_graph_type: str = '00',
    n_jobs=1
):

    multiplex, bipartite, eta, lamb = init_GRN(
        TF_layer=TF_layer,
        ATAC_layer=ATAC_layer,
        RNA_layer=RNA_layer,
        TF_ATAC_bipartite=TF_ATAC_bipartite,
        ATAC_RNA_bipartite=ATAC_RNA_bipartite,
        TF_layer_graph_type=TF_layer_graph_type,
        ATAC_layer_graph_type=ATAC_layer_graph_type,
        RNA_layer_graph_type=RNA_layer_graph_type,
        TF_ATAC_bipartite_graph_type=TF_ATAC_bipartite_graph_type,
        ATAC_RNA_bipartite_graph_type=ATAC_RNA_bipartite_graph_type
        )
    print("initialization done.")

    if type(seeds) is str:
        seeds = [seeds]
        print(seeds)

    ranking = hummuspy.explore_network.compute_multiple_RandomWalk(
        multiplex,
        bipartite,
        eta=eta.tolist(),
        lamb=numpy.array(lamb).tolist(),
        seeds=seeds,
        save=False,
        spec_layer_result_saved='RNA',
        n_jobs=n_jobs)

    return ranking


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
    df : pandas.DataFrame
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

    # process the config to the right format to
    # compute the random walks without local saving
    config = hummuspy.config.process_config(config, multilayer_f)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            gene_list = numpy.unique(numpy.concatenate([gene_list, layer_nodes]
                                                       ))

    config['seeds'] = gene_list

    df = hummuspy.explore_network.compute_multiple_RandomWalk(
        **config,
        output_f=output_f,
        save=False,
        return_df=return_df,
        spec_layer_result_saved=tf_multiplex,
        n_jobs=njobs)

    df['gene'] = df['seed']
    df['tf'] = df['target']
    del df['target']
    del df['seed']

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            tf_list = numpy.unique(numpy.concatenate([tf_list, layer_nodes]))
        tf_list = tf_list[tf_list != 'fake_node']

    # Add normalisation ?
    df = df[df['tf'].isin(tf_list)]

    if save is True:
        print('Saving the result in ', output_f)
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
    df : pandas.DataFrame
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

    # process the config to the right format to
    # compute the random walks without local saving
    config = hummuspy.config.process_config(config, multilayer_f)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            gene_list = numpy.unique(numpy.concatenate([gene_list, layer_nodes]
                                                       ))
    config['seeds'] = gene_list

    df = hummuspy.explore_network.compute_multiple_RandomWalk(
        **config,
        output_f=output_f,
        save=False,
        return_df=return_df,
        spec_layer_result_saved=peak_multiplex,
        n_jobs=njobs)

    df['gene'] = df['seed']
    df['peak'] = df['target']
    del df['target']
    del df['seed']

    if peak_list is None:
        peak_list = []
        for layer in config['multiplex'][peak_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            peak_list = numpy.unique(numpy.concatenate([peak_list, layer_nodes]
                                                       ))

    # Add normalisation ?
    df = df[df['peak'].isin(peak_list)]

    if save is True:
        print('Saving the result in ', output_f)
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
    df : pandas.DataFrame
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

    # process the config to the right format to
    # compute the random walks without local saving
    config = hummuspy.config.process_config(config, multilayer_f)

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            tf_list = numpy.unique(numpy.concatenate([tf_list, layer_nodes]
                                                     ))
        tf_list = tf_list[tf_list != 'fake_node']

    config['seeds'] = tf_list

    df = hummuspy.explore_network.compute_multiple_RandomWalk(
        **config,
        output_f=output_f,
        save=False,
        return_df=return_df,
        spec_layer_result_saved=peak_multiplex,
        n_jobs=njobs)

    df['tf'] = df['seed']
    df['peak'] = df['target']
    del df['target']
    del df['seed']

    if peak_list is None:
        peak_list = []
        for layer in config['multiplex'][peak_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            peak_list = numpy.unique(numpy.concatenate([peak_list, layer_nodes]
                                                       ))

    # Add normalisation ?
    df = df[df['peak'].isin(peak_list)]

    if save is True:
        print('Saving the result in ', output_f)
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the binding regions prediction result.'
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
    df : pandas.DataFrame
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

    # process the config to the right format to
    # compute the random walks without local saving
    config = hummuspy.config.process_config(config, multilayer_f)

    if gene_list is None:
        gene_list = []
        for layer in config['multiplex'][rna_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            gene_list = numpy.unique(numpy.concatenate([gene_list, layer_nodes]
                                                       ))

    if tf_list is None:
        tf_list = []
        for layer in config['multiplex'][tf_multiplex]['layers']:
            df_layer = pandas.read_csv(
                layer,
                sep='\t',
                header=None,
                index_col=None)

            layer_nodes = numpy.concatenate([
                numpy.unique(df_layer[0].values),
                numpy.unique(df_layer[1].values)])
            tf_list = numpy.unique(numpy.concatenate([tf_list, layer_nodes]
                                                     ))
        tf_list = tf_list[tf_list != 'fake_node']

    config['seeds'] = tf_list

    df = hummuspy.explore_network.compute_multiple_RandomWalk(
        **config,
        output_f=output_f,
        save=False,
        return_df=return_df,
        spec_layer_result_saved=rna_multiplex,
        n_jobs=njobs)

    df['tf'] = df['seed']
    df['gene'] = df['target']
    del df['target']
    del df['seed']

    # Add normalisation ?
    df = df[df['gene'].isin(gene_list)]

    if save is True:
        print('Saving the result in ', output_f)
        assert output_f is not None, 'You need to provide an output_f name ' +\
            'to save the terget gene predictions result'
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
    df : pandas.DataFrame
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
