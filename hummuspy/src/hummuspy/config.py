import typing

import yaml

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
    
def group_per_layer(
        multiplex_list
        ):
    """Group multiplex info per layer.
    
    Structure returned:
    {multiplex_name: {'layers': [layer1, layer2, ...],
                      'graph_type': [graph_type1, graph_type2, ...]}
     ...}
    
    """
    multiplex_organised = dict()
    for multiplex_name in list(multiplex_list['networks'].keys()):
        # we instanciata a dict for each multiplex network
        multiplex_organised[multiplex_name] = dict()
        # we add the layer names
        multiplex_organised[multiplex_name]['layers'] =\
            [layer for layer in multiplex_list['networks'][network_name]]
        # we add the graph type
        multiplex_organised[multiplex_name]['graph_type'] =\
            [str(int(directed))+str(int(weighted)) for directed, weighted in\
             zip(multiplex_list['directed'][multiplex_name],\
                 multiplex_list['weighted'][multiplex_name])]
    return multiplex_organised




def general_config(
    multiplexes: typing.Union[dict[dict[str]]],
    bipartites: typing.Union[str, list[str], dict[str]],
    seed_path: str = 'seeds/seeds.txt',
    folder_multiplexes = 'multiplex',
    folder_bipartites = 'bipartite',
    bipartite_type: typing.Union[str, list[str], dict[str]] = 'undirected',
    self_loops = 0,
    restart_prob = 0.7,):
    
    """Create a very general config file for the hummus pipeline."""
    config = dict()
    config['multiplex'] = dict()
    config['bipartite'] = dict()
    config['seed'] = seed_path
    config['self_loops'] = self_loops

    # We add the multiplexes to the config
    for multiplex_name in multiplexes:
        # If folder_multiplexes is None we use the multiplex name as folder name
        config['multiplex'][multiplex_name] = dict()
        config['multiplex'][multiplex_name]['layers'] =\
            [(folder_multiplexes+'/'+multiplex_name+'/'+layer).replace('//', '/') for layer in multiplexes[multiplex_name]['layers']]
        config['multiplex'][multiplex_name]['graph_type'] =\
            multiplexes[multiplex_name]['graph_type']
        

    # if type of bipartites not associated to their names already,
    # we create a dict with the same order as the bipartites
    bipartite_type = make_values_list(bipartite_type)
    if type(bipartite_type) == list:
        temp = dict()
        for i in range(len(bipartites)):
            temp[bipartites.keys()[i]] = bipartite_type[i]
        bipartite_type = temp

    # we add the bipartites
    for bipartite in bipartites:
        bipartite_loc = folder_bipartites+'/'+bipartite
        config['bipartite'][bipartite_loc] = dict()
        config['bipartite'][bipartite_loc]['source'] = bipartites[bipartite]['multiplex_left']
        config['bipartite'][bipartite_loc]['target'] = bipartites[bipartite]['multiplex_right']
        config['bipartite'][bipartite_loc]['graph_type'] = bipartite_type[bipartite]

    config['r'] = restart_prob
    return config

def save_config(config, filename):
    with open(filename, 'w') as f:
        yaml.dump(config, f)

