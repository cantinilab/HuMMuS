from typing import Union
import yaml
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx


def make_values_list(
        values,
        types: Union[str, bool, int, float] = (str, bool, int, float)):
    """Transform layer type to be sure a list of values is returned.
    Parameters
    ----------
    values: list, dict, str, bool, int, float
        The values to transform into a list.
    types: list, dict, str, bool, int, float
        The types of the values to transform into a list.
    Returns
    -------
    values: list
        The values transformed into a list.
    Raises
    ------
    TypeError if the values are not of the given types.

    """
    types = list(types)
    if type(values) == list:
        return values
    elif type(values) == tuple:
        return list(values)
    elif type(values) == dict:
        return list(values.values())
    elif type(values) in types:
        return [values]
    else:
        raise TypeError('Layer name(s) {} should be given through {},'.format(
            values,
            types) +
            'a list, or a dictionary(where the key would be the file_paths)')


def group_per_layer(
        multiplex_list
        ):
    """Group multiplex info per layer in a dictionary,
    from a dictionary of multiplexes obtained from hummus R objects.

    Parameters
    ----------
    multiplex_list: dict
        Dictionary of multiplexes obtained from hummus R objects.
        Structure:
        {multiplex1: {layer1: graph_type1, layer2: graph_type2, ...},
        multiplex2: {layer1: graph_type1, layer2: graph_type2, ...},
        ...}

    Returns
    -------
    multiplex_organised: dict
        Dictionary of multiplexes organised in 2 lists,
        one for the layer names, one for the graph types.
        Structure:
        {multiplex1: {layers: [layer1, layer2, ...],
                      graph_type: [graph_type1, graph_type2, ...]},
        multiplex2: {layers: [layer1, layer2, ...],
                     graph_type: [graph_type1, graph_type2, ...]},
        ...}
    """
    print(multiplex_list)
    multiplex_organised = dict()
    for multiplex_name in multiplex_list:
        # we instanciate a dict for each multiplex network
        multiplex_organised[multiplex_name] = dict()
        # we add the layer names
        multiplex_organised[multiplex_name]['layers'] =\
            [layer for layer in multiplex_list[multiplex_name]]
        # we add the graph type
        multiplex_organised[multiplex_name]['graph_type'] =\
            [multiplex_list[multiplex_name][layer]
                for layer in multiplex_list[multiplex_name]]
    return multiplex_organised


def general_config(
        multiplexes: dict[dict[str]],
        bipartites: Union[str, list[str], dict[str]],
        seed_path: str = 'seeds/seeds.txt',
        folder_multiplexes='multiplex',
        folder_bipartites='bipartite',
        self_loops=0,
        restart_prob=0.7,
        bipartites_type: Union[str, list[str], dict[str]] = ('00', '00'),
        save_configfile: bool = True,
        config_filename: str = 'config.yaml',
        suffix='.tsv'):

    """Create a very general config file for the hummus pipeline.
    The config file is a dictionary that can be saved as a yaml file.
    The config file is organised as follows:
    config = {'multiplex': {'multiplex1': {'layers': [layer1, layer2, ...],
                                           'graph_type': [graph_type1, ...]},
                            'multiplex2': {'layers': [layer1, layer2, ...],
                                           'graph_type': [graph_type1, ...]},
                            ...},
              'bipartite': {'bipartite1': {'source': multiplex1,
                                           'target': multiplex2},
                            'bipartite2': {'source': multiplex1,
                                           'target': multiplex2},
                                           ...},
               'seed': seed_path,
               'self_loops': self_loops,
               'r': restart_prob}

    Parameters
    ----------
    multiplexes: dict
        Dictionary of multiplexes obtained from hummus R objects.
        Structure:
        {multiplex1: {layer1: graph_type1, layer2: graph_type2, ...},
        multiplex2: {layer1: graph_type1, layer2: graph_type2, ...},
        ...}
    bipartites: str, list[str], dict[str]
        The names of the bipartites.
        If str, the name of the bipartite.
        If list[str], the names of the bipartites.
        If dict[str], the names of the bipartites as keys and
        the type of the bipartites as values.
        e.g.: bipartites = {'bipartite1': {multiplex_right: multiplex1,
                                           multiplex_left: multiplex2},
                            'bipartite2': {multiplex_right: multiplex1,
                                           multiplex_left: multiplex2}}
    seed_path: str
        The path to the seed file.
    folder_multiplexes: str
        The folder where the multiplexes are located.
    folder_bipartites: str
        The folder where the bipartites are located.
    self_loops: int
        The number of self loops to add to the multiplexes.
    restart_prob: float
        The restart probability for the RWR.
    bipartites_type: str, list[str], dict[str]
        The type of the bipartites.
        If str, the type of the bipartite.
        If list[str], the types of the bipartites.
        If dict[str], the names of the bipartites as keys and
        the type of the bipartites as values.
        e.g.: bipartites_type = {'bipartite1': '00',
                                 'bipartite2': '00'}
    save_configfile: bool
        If True, save the config file as a yaml file.
    config_filename: str
        The name of the config file to save.
    suffix: str
        The suffix of the multiplex and bipartite files.

    Returns
    -------
    config: dict
        The config dictionary.

    e.g.:
    config = {
            'multiplex': {
                'RNA': {'layers': ['multiplex/RNA/RNA_GENIE3.tsv'],
                        'graph_type': ['00']},
                'TF': {'layers': ['multiplex/TF/TF_network.tsv'],
                       'graph_type': ['00']},
                'peaks': {'layers': ['multiplex/peaks/peak_network_G3.tsv'],
                          'graph_type': ['00']}},
            'bipartite': {'bipartite/RNA_peak': {'source': 'peaks',
                                                 'target': 'TF',
                                                 'graph_type': '00'},
                          'bipartite/TF_peak': {'source': 'peaks',
                                                'target': 'RNA',
                                                'graph_type': '00'}},
            'seed': 'seeds/seeds.txt',
            'self_loops': 0,
            'r': 0.7,
            'eta': [1, 0, 0],
            'lambda':  [[0.5, 0.0, 0.5],
                        [0.0, 0.5, 0.5],
                        [1/3, 1/3, 1/3]]}

    Raises
    ------
    AssertionError if the length of eta is not equal to the number of layers.
    """

    multiplexes = group_per_layer(multiplexes)

    config = dict()
    config['multiplex'] = dict()
    config['bipartite'] = dict()
    config['seed'] = seed_path
    config['self_loops'] = self_loops

    # We add the multiplexes to the config
    for multiplex_name in multiplexes:
        # If folder_multiplexes is None, use the multiplex name as folder name
        config['multiplex'][multiplex_name] = dict()
        config['multiplex'][multiplex_name]['layers'] =\
            [(folder_multiplexes+'/'+multiplex_name+'/'+layer+suffix)
             .replace('//', '/')
             for layer in multiplexes[multiplex_name]['layers']]
        config['multiplex'][multiplex_name]['graph_type'] =\
            multiplexes[multiplex_name]['graph_type']

    # if type of bipartites not associated to their names already,
    # we create a dict with the same order as the bipartites
    if type(bipartites_type) == list or type(bipartites_type) == tuple:
        print("bipartites_type has been provided throguh a list, make sure " +
              "the order matches the one of the 'bipartites' dictionary' keys."
              )
        temp = dict()
        for i in range(len(bipartites)):
            temp[list(bipartites.keys())[i]] = bipartites_type[i]
        bipartites_type = temp
    elif type(bipartites_type) == dict:
        assert list(bipartites.keys()).sort() == list(
            bipartites_type.keys()).sort(),\
            "The keys of the 'bipartites_type' and of the 'bipartites' " +\
            "dictionary doesn't seem to match. " + \
            "Please provide identical keys for each dictionary."
    else:
        raise TypeError("bipartites_type should be a list, a tuple or a dictionary")

    # we add the bipartites
    for bipartite in bipartites:
        bipartite_loc = folder_bipartites+'/'+bipartite
        config['bipartite'][bipartite_loc] = dict()
        config['bipartite'][bipartite_loc]['source'] = \
            bipartites[bipartite]['multiplex_left']
        config['bipartite'][bipartite_loc]['target'] = \
            bipartites[bipartite]['multiplex_right']
        config['bipartite'][bipartite_loc]['graph_type'] = \
            bipartites_type[bipartite]
    config['r'] = restart_prob

    if save_configfile is True:
        save_config(config, config_filename)

    return config


def save_config(config, filename):
    """ Save the config dictionary as a yaml file.
    The parser used is the default one of yaml python package.

    Parameters
    ----------
    config: dict
        The config dictionary.
    filename: str
        The name of the file to save the config to.

    Returns
    -------
    None
    """

    config = dict(config)
    multiplex_order = list(config['multiplex'].keys())
    # we sort eta, lamb and multiplexes to make sure the order is always the same
    multiplex_order.sort()
    config['multiplex'] = {k: config['multiplex'][k] for k in multiplex_order}
    if 'eta' in config.keys():
        config['eta'] = config['eta'].loc[multiplex_order]
    if 'lamb' in config.keys():
        config['lamb'] = config['lamb'].loc[multiplex_order, multiplex_order]

    config = setup_proba_config(config, lamb=config["lamb"], eta=config["eta"])

    with open(filename, 'w') as f:
        yaml.dump(config, f)


def open_config(filename):
    with open(filename) as file:
        config_dic = yaml.load(file, Loader=yaml.BaseLoader)

    if "lamb" in config_dic.keys():
        config_dic["lamb"] = pd.DataFrame(
            config_dic["lamb"],
            index=list(config_dic["multiplex"].keys()),
            columns=list(config_dic["multiplex"].keys())
            ).astype(float)
    if "eta" in config_dic.keys():
        config_dic["eta"] = pd.Series(
            config_dic["eta"],
            index=list(config_dic["multiplex"].keys())
            ).astype(float)

    return config_dic


def old_setup_proba_config(
        config: dict,
        eta: Union[list[float], pd.Series, np.ndarray],
        lamb: pd.DataFrame):
    """ Setup the RWR probability for the exploration of hummus networks
    with the given eta and lambda values.
    The lambda values are normalised (per rows) to sum to 1.

    Parameters
    ----------
    config: dict
        The config dictionary.
    eta: list[float]
        The eta values for the RWR probability, must sum up to 1.
        e.g.: [0.5, 0.5, 0]

    lamb: list[list[float]]
        The lambda values for the RWR probability.
        e.g.: lamb = pd.DataFrame(np.ones((3,3)),
                           index = ['TF', 'Gene', 'Peak''],
                           columns = ['TF', 'Gene', 'Peak'])
        lamb.loc[i, j] corresponds to the probability
        to go from layer i to layer j.

    Returns
    -------
    config: dict
        The config dictionary with the RWR probability setup.

    Raises
    ------
    AssertionError if the length of eta is not equal to the number of layers.

    e.g.:
    config = {'multiplex': {'TF': {'file_path': 'TF.csv',
                                   'type': 'TF'},
                            'Gene': {'file_path': 'Gene.csv',
                                     'type': 'Gene'},
                            'Peak': {'file_path': 'Peak.csv',
                                     'type': 'Peak'}},
              'bipartite': {'TF_Peak': {'source': 'TF',
                                        'target': 'Peak'},
                            'Gene_Peak': {'source': 'Gene',
                                          'target': 'Peak'}},
              'eta': [0.5, 0.5, 0],
              'lambda': [[0.5, 0.5, 0],
                         [0.5, 0, 0.5],
                         [0, 0.5, 0.5]]}
    """

    # Check that the length of eta is equal to the number of layers
    assert len(config['multiplex']) == len(eta),\
        "eta (length of {}) should be as long as the number of layers ({})"\
        .format(len(eta), len(config['multiplex']))

    # Normalise lamb per rows
    lamb = lamb.div(lamb.sum(axis=1), axis=0)
    # Check that lamb is a valid probability matrix
    assert check_lamb(lamb, config),\
        "lamb is not a valid probability matrix according to bipartites"
    # Transform eta to list if it's a pandas Series or a numpy array
    if type(eta) == pd.Series:
        eta = eta.values.tolist()
    elif type(eta) == np.ndarray:
        eta = eta.tolist()

    # Add eta and lamb to config
    config['eta'] = eta
    config['lamb'] = lamb.values.tolist()

    return config


def setup_proba_config(
        config: dict,
        eta: Union[pd.Series, np.ndarray, list[float]],
        lamb: Union[pd.DataFrame, np.ndarray]
        ):
    """ Setup the RWR probability for the exploration of hummus networks
    with the given eta and lambda values.
    The lambda values are normalised (per columns) to sum to 1.

    Parameters
    ----------
    config: dict
        The config dictionary.
    eta: Union[pd.Series, np.array]
        The eta values for the RWR probability, must sum up to 1.
        e.g.: pd.Series([0.5, 0.5, 0], index = ['TF', 'peaks', 'RNA'])
        If a pandas Series is provided, the index must be the layer names.
        If a numpy array is provided, the order of the values must be
        the same as the order of the layers in the config.

    lamb: pd.DataFrame
        The lambda values for the RWR probability.
        e.g.: lamb = pd.DataFrame(np.ones((3,3)),
                           index = ['TF', 'Gene', 'Peak''],
                           columns = ['TF', 'Gene', 'Peak'])
        lamb.loc[i, j] corresponds to the probability
        to go from layer j to layer i.

    Returns
    -------
    config: dict
        The config dictionary with the RWR probability setup.

    Raises
    ------
    AssertionError if the length of eta is not equal to the number of layers.

    e.g.:
    config = {'multiplex': {'TF': {'file_path': 'TF.csv',
                                   'type': 'TF'},
                            'Gene': {'file_path': 'Gene.csv',
                                     'type': 'Gene'},
                            'Peak': {'file_path': 'Peak.csv',
                                     'type': 'Peak'}},
              'bipartite': {'TF_Peak': {'source': 'TF',
                                        'target': 'Peak'},
                            'Gene_Peak': {'source': 'Gene',
                                          'target': 'Peak'}},
              'eta': [0.5, 0.5, 0],
              'lambda': [[0.5, 0.5, 0],
                         [0.5, 0, 0.5],
                         [0, 0.5, 0.5]]}
    """
    if type(eta) == np.array or type(eta) == list:
        # Check that the length of eta is equal to the number of layers
        assert len(config['multiplex']) == len(eta),\
            "eta (length : {}) should be as long as the number of layers ({})"\
            .format(len(eta), len(config['multiplex']))
        # Transform eta to pandas Series
        eta = pd.Series(eta, index=config['multiplex'].keys())
    elif type(eta) == pd.Series:
        # Check that the index of eta are the layer names
        assert eta.index.tolist().sort() == list(
            config['multiplex'].keys()).sort(),\
            "eta index ({}) should be the same as the layer names ({})"\
            .format(eta.index.tolist(), list(config['multiplex'].keys()))
    else:
        raise TypeError("eta should be a numpy array or a pandas Series")

    if type(lamb) == np.array:
        # Check that lamb is a square matrix of size len(config['multiplex'])
        assert lamb.shape[0] == lamb.shape[1],\
            "lamb should be a square matrix"
        assert lamb.shape[0] == len(config['multiplex']),\
            "lamb should be a square matrix of size {}".format(
            len(config['multiplex']))
        # Transform lamb to pandas DataFrame
        lamb = pd.DataFrame(lamb, index=config['multiplex'].keys(),
                            columns=config['multiplex'].keys())
    elif type(lamb) == pd.DataFrame:
        # Check that lamb index and columns are the layer names
        assert lamb.index.tolist().sort() == list(
            config['multiplex'].keys()).sort(),\
            "lamb index ({}) should be the same as the layer names ({})"\
            .format(lamb.index.tolist(), list(config['multiplex'].keys()))
        assert lamb.columns.tolist().sort() == list(
            config['multiplex'].keys()).sort(),\
            "lamb columns ({}) should be the same as the layer names ({})"\
            .format(lamb.columns.tolist(), list(config['multiplex'].keys()))
    else:
        raise TypeError("lamb should be a numpy array or a pandas DataFrame")

#    lamb = lamb.transpose()
    # Normalise lamb per col
    lamb = lamb.div(lamb.sum(axis=0), axis=1)
    lamb = lamb.fillna(0)

    # Check that lamb is a valid probability matrix
    assert check_lamb(lamb, config),\
        "lamb is not a valid probability matrix according to bipartites"

    # Order eta and lamb according to the order of the layers in the config
    eta = eta[config['multiplex'].keys()]
    lamb = lamb.loc[config['multiplex'].keys(), config['multiplex'].keys()]

    # Transform eta to list if it's a pandas Series or a numpy array
    if type(eta) == pd.Series:
        eta = eta.values.tolist()
    elif type(eta) == np.ndarray:
        eta = eta.tolist()

    # Add eta and lamb to config
    config['eta'] = eta
    config['lamb'] = lamb.values.tolist()

    return config


def initialise_lamb(config,
                    tf_multiplex='TF',
                    peak_multiplex='peaks',
                    rna_multiplex='RNA',
                    value=1):

    for multiplex in [tf_multiplex, peak_multiplex, rna_multiplex]:
        assert multiplex in config['multiplex'].keys(),\
            "The multiplex {}".format(multiplex) +\
            " is not in the config file provided"

    ordered_multiplexes = config['multiplex'].keys()

    if value == 1:
        array = np.ones((len(ordered_multiplexes), len(ordered_multiplexes)))
    elif value == 0:
        array = np.zeros((len(ordered_multiplexes), len(ordered_multiplexes)))
    else:
        raise ValueError('value param of initialise_lamb should be 1 or 0')

    lamb = pd.DataFrame(array,
                        index=ordered_multiplexes,
                        columns=ordered_multiplexes)
    return lamb


def get_single_layer_eta(config, starting_multiplex='RNA'):

    ordered_multiplex = config['multiplex'].keys()
    assert starting_multiplex in ordered_multiplex,\
        "It seems starting_multiplex not in config['multiplex']"

    eta = pd.Series([0 for k in ordered_multiplex],
                    index=ordered_multiplex)
    eta[starting_multiplex] = 1
    return eta


#############################
# 1/4 Get GRN classic lamb  #
#############################
def get_grn_lamb(config,
                 tf_multiplex='TF',
                 peak_multiplex='peaks',
                 rna_multiplex='RNA',
                 draw=False
                 ):

    lamb = initialise_lamb(config,
                           tf_multiplex,
                           peak_multiplex,
                           rna_multiplex,
                           value=1)  # because enhancer lamb is mostly 1s

    # Remove proba between TF and RNA layers
    lamb.loc[tf_multiplex, rna_multiplex] = 0
    lamb.loc[rna_multiplex, tf_multiplex] = 0
    lamb = lamb.transpose()
    lamb = lamb.div(lamb.sum(axis=0),
                    axis=1)

    # max_lambd check to see
    assert check_lamb(config=config, lamb=lamb), "There seem to be a " +\
        "incoherence between bipartite source/targets and multiplex names" +\
        "provided in get_classic_grn_lamb"

    if draw is True:
        to_draw_lamb = lamb.loc[[tf_multiplex, peak_multiplex, rna_multiplex],
                                [tf_multiplex, peak_multiplex, rna_multiplex]]

        draw_lamb(to_draw_lamb)

    return lamb


##################################
# 2/4 Get enhancers classic lamb #
##################################
def get_enhancers_lamb(config,
                       tf_multiplex='TF',
                       peak_multiplex='peaks',
                       rna_multiplex='RNA',
                       draw=False
                       ):

    lamb = initialise_lamb(config,
                           tf_multiplex,
                           peak_multiplex,
                           rna_multiplex,
                           value=0)  # because enhancer lamb is mostly 0s

    # Add proba between peaks and RNA layers
    lamb.loc[peak_multiplex, peak_multiplex] = 1
    lamb.loc[rna_multiplex, peak_multiplex] = 1
    lamb.loc[tf_multiplex, peak_multiplex] = 1

    lamb = lamb.transpose()
    lamb = lamb.div(lamb.sum(axis=0),
                    axis=1)
    lamb = lamb.fillna(0)

    # max_lambd check to see
    assert check_lamb(config=config, lamb=lamb), "There seem to be a " +\
        "incoherence between bipartite source/targets and multiplex names" +\
        "provided in get_classic_grn_lamb"

    if draw is True:
        to_draw_lamb = lamb.loc[[tf_multiplex, peak_multiplex, rna_multiplex],
                                [tf_multiplex, peak_multiplex, rna_multiplex]]
        draw_lamb(to_draw_lamb)

    return lamb


########################################
# 3/4 Get binding regions classic lamb #
########################################
def get_binding_regions_lamb(config,
                             tf_multiplex='TF',
                             peak_multiplex='peaks',
                             rna_multiplex='RNA',
                             draw=False
                             ):

    lamb = initialise_lamb(config,
                           tf_multiplex,
                           peak_multiplex,
                           rna_multiplex,
                           value=0)  # because enhancer lamb is mostly 0s

    # Add proba between TF and peaks layers
    lamb.loc[tf_multiplex, tf_multiplex] = 1
    lamb.loc[tf_multiplex, peak_multiplex] = 1
    lamb.loc[peak_multiplex, peak_multiplex] = 1
    lamb.loc[peak_multiplex, tf_multiplex] = 1
    lamb.loc[rna_multiplex, peak_multiplex] = 1

    lamb = lamb.transpose()
    lamb = lamb.div(lamb.sum(axis=0),
                    axis=1)
    lamb = lamb.fillna(0)

    # max_lambd check to see
    assert check_lamb(config=config, lamb=lamb), "There seem to be a " +\
        "incoherence between bipartite source/targets and multiplex names" +\
        "provided in get_classic_grn_lamb"

    if draw is True:
        to_draw_lamb = lamb.loc[[tf_multiplex, peak_multiplex, rna_multiplex],
                                [tf_multiplex, peak_multiplex, rna_multiplex]]
        draw_lamb(to_draw_lamb)

    return lamb


########################################
# 4/4 Get binding regions classic lamb #
########################################
def get_target_genes_lamb(config,
                          tf_multiplex='TF',
                          peak_multiplex='peaks',
                          rna_multiplex='RNA',
                          draw=False
                          ):

    lamb = initialise_lamb(config,
                           tf_multiplex,
                           peak_multiplex,
                           rna_multiplex,
                           value=1)  # because enhancer lamb is mostly 0s

    # Remove proba between TF and RNA layers
    lamb.loc[tf_multiplex, rna_multiplex] = 0
    lamb.loc[rna_multiplex, tf_multiplex] = 0
    lamb.loc[peak_multiplex, tf_multiplex] = 0  # can't go back up to TF

    lamb = lamb.transpose()
    lamb = lamb.div(lamb.sum(axis=0),
                    axis=1)
    lamb = lamb.fillna(0)

    # max_lambd check to see
    assert check_lamb(config=config, lamb=lamb), "There seem to be a " +\
        "incoherence between bipartite source/targets and multiplex names" +\
        "provided in get_classic_grn_lamb"

    if draw is True:
        to_draw_lamb = lamb.loc[[tf_multiplex, peak_multiplex, rna_multiplex],
                                [tf_multiplex, peak_multiplex, rna_multiplex]]
        draw_lamb(to_draw_lamb)

    return lamb


##############################
# Check proba transitions    #
##############################
def get_max_lamb(config, directed=True, draw=False, figsize=(7, 7)):
    """Calculate the maximum lamb matrix according to bipartites
    indicated in the config in input.
    Bipartites are used to know which layers can be connected to each other.

    Parameters
    ----------
    config: dict
        The config dictionary.

    Returns
    -------
    max_lamb: pandas.DataFrame
        The maximum lamb matrix according to bipartites.
        Structure:
          pd.DataFrame([[0.5, 0.5, 0],
                        [1/3, 1/3, 1/3],
                        [0, 0.5, 0.5]],
                        index = ['TF', 'Peak', 'RNA'],
                        columns = ['TF', 'Peak', 'RNA'])
    """
    # Get layer names connected by bipartites

    bipartites = dict(config['bipartite'])
    bip_keys = list(bipartites.keys())
    for bipartite in bip_keys:
        # if the bipartite is undirected, we add the inversed bipartite
        if bipartites[bipartite]["graph_type"][0]=="0":
            bipartites[bipartite+"_inversed"] = {
                "source": bipartites[bipartite]["target"],
                "target": bipartites[bipartite]["source"],
                "graph_type": bipartites[bipartite]["graph_type"]
                }
            
    if directed != 'inversed':
        positions_bipartites = [(bipartites[bipartite]['target'],
                                 bipartites[bipartite]['source'])
                                 for bipartite in bipartites]
    else: # if we want to inverse source and target layers
        positions_bipartites = [(bipartites[bipartite]['source'],
                                 bipartites[bipartite]['target'])
                                 for bipartite in bipartites]

    # Create an empty dataframe with layer names as index and columns
    # that we'll fill where it's possible according to the bipartites
    max_lamb = pd.DataFrame(np.zeros((len(config['multiplex']),
                                      len(config['multiplex']))),
                            index=config['multiplex'].keys(),
                            columns=config['multiplex'].keys())

    # fill each position corresponding to bipartite options
    for position in positions_bipartites:
        max_lamb.loc[position] = 1
    # Future : Since we can inverse source and target layers)
    # could be conditionned by bipartite directionality
    if directed:
        pass
    else:
        max_lamb += max_lamb.transpose()
    # filling the diagonal to allow intra-layer exploration
    max_lamb += np.eye(len(config['multiplex'])).astype(int)
    # normalise per rows
    max_lamb = max_lamb.div(max_lamb.sum(axis=0), axis=1)

    if draw is True:
        draw_lamb(max_lamb, figsize=figsize)

    return max_lamb


def check_lamb(lamb, config, directed=True, draw=False):
    """Check that lamb is a valid probability matrix according to bipartites
    indicated in the config in input.
    Bipartites are used to know which layers can be connected to each other.

    Parameters
    ----------
    lamb: pandas.DataFrame
        The lamb matrix.
    config: dict
        The config dictionary.

    Returns
    -------
    bool
        True if lamb is a valid probability matrix according to bipartites.
        False otherwise.
    """

    # Calculate the maximum lamb matrix according to bipartites
    max_lamb = (get_max_lamb(config, directed=directed, draw=draw) > 0).astype(int)

    # Count number of locations where lamb has no non-zero value
    # and max_lamb is zero
    X = np.sum(np.sum((max_lamb - lamb) < 0, axis=0), axis=0)

    # If X > 0, lamb is not a valid probability matrix according to bipartites.
    if X > 0:
        return False
    # Else, lamb is a valid probability matrix according to bipartites.
    else:
        return True
    


##############################
# Draw proba transitions     #
##############################
def draw_config(config, figsize=(7, 7)):
    """Draw the lamb argument config as a networkx graph.
    
    Parameters
    ----------
    config: dict
        The config dictionary.
    figsize: A tuple containing the dimensions of the plot.

    Returns
    -------
    None
    """

    # Base node color on starting multiplex (from eta)
    node_color = config["eta"][config['lamb'].columns].values
    node_color = ["blue" if color == 0 else "green"
                  for color in node_color]

    draw_lamb(config["lamb"], figsize=figsize, node_color=node_color)


def draw_lamb(df, figsize=(7, 7), node_color='blue'):
    """Draw the lamb matrix as a networkx graph.

    Parameters
    ----------
    df: pandas.DataFrame
        The lamb matrix.
        Structure:
            pd.DataFrame([[0.5, 0.5, 0],
                            [1/3, 1/3, 1/3],
                            [0, 0.5, 0.5]],
                            index = ['TF', 'Peak', 'RNA'],
                            columns = ['TF', 'Peak', 'RNA'])

    figsize: A tuple containing the dimensions of the plot.
    node_color: Union[str, list[str]]
        The color of the nodes.
        If str, the color of all the nodes.
        If list[str], colors of each node.
        e.g.: node_color = 'blue'
              node_color = ['blue', 'red', 'blue']

    Returns
    -------
    None

    Examples
    --------
    >>> lamb = pd.DataFrame([[0.5, 0.5, 0],
                                [1/3, 1/3, 1/3],
                                [0, 0.5, 0.5]],
                                index = ['TF', 'Peak', 'RNA'],
                                columns = ['TF', 'Peak', 'RNA'])
    >>> draw_lamb(lamb)

    See Also
    --------
    draw_networkx
    """

    # Create a figure and an axis
    fig, ax = plt.subplots(figsize=figsize)

    # Create a directed graph from the lamb matrix
    G = nx.from_pandas_adjacency(df.transpose(), create_using=nx.DiGraph())
    pos = {list(G.nodes)[-i-1]: np.array((0, i))
           for i in range(-len(G.nodes), 0)}

    # Draw the graph without edge labels
    nx.draw_networkx(G,
                     with_labels=True,
                     pos=pos,
                     node_size=1500,
                     width=4,
                     node_color=node_color,
                     alpha=0.8,
                     font_weight="bold",
                     arrows=True,
                     connectionstyle='arc3, rad = 0.6')

    # Get edge weights that will be used as edge labels
    edge_labels = nx.get_edge_attributes(G, 'weight')
    # Round edge weights to 2 decimals
    edge_labels = {k: round(v, 2) for k, v in edge_labels.items()}

    # Separate self edges from non self edges
    self_edges = {}
    non_self_edges = {}
    for edge in edge_labels.keys():
        # If the edge is a self edge we add it to self_edges
        if edge[0] == edge[1]:
            self_edges[edge] = edge_labels[edge]
        # If the edge is not a self edge we add it to non_self_edges
        else:
            non_self_edges[edge] = edge_labels[edge]

    # Draw self edges labels
    # We add 0.3 to the y coordinate to avoid overlapping with the node
    pos_self_labels = {k: np.array([0, pos[k][1]+0.30]) for k in pos}
    my_draw_networkx_edge_labels(G,
                                 pos_self_labels,
                                 edge_labels=self_edges,
                                 font_color='k',
                                 font_size=12,
                                 label_pos=12,
                                 rad=0.6,
                                 rotate=False)

    # Draw non self edges labels
    my_draw_networkx_edge_labels(G,
                                 pos,
                                 edge_labels=non_self_edges,
                                 font_color='k',
                                 font_size=12,
                                 label_pos=0,
                                 rad=0.6,
                                 rotate=False)

    ax.set_ylim([list(pos.values())[0][1] - 0.5,
                 list(pos.values())[-1][1] + 0.5])

def my_draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0
):
    """Draw edge labels.

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center',
                            'top',
                            'bottom',
                            'baseline',
                            'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5*pos_1 + 0.5*pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0, 1),
                                    (-1, 0)])
        ctrl_1 = linear_mid + rad*rotation_matrix@d_pos
        ctrl_mid_1 = 0.5*pos_1 + 0.5*ctrl_1
        ctrl_mid_2 = 0.5*pos_2 + 0.5*ctrl_1
        bezier_mid = 0.5*ctrl_mid_1 + 0.5*ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round",
                        ec=(1.0, 1.0, 1.0),
                        fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items
