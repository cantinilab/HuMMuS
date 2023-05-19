import typing
import yaml

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import networkx as nx


def make_values_list(
        values,
        types: typing.Union[str, bool, int, float]):
    """Transform layer type to be sure a list of values is returned."""
    if type(values) == list:
        return values
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

    with open(filename, 'w') as f:
        yaml.dump(config, f)


def setup_proba_config(
        config: dict,
        eta: typing.Union[list[float], pd.Series, np.ndarray],
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
    config['lambda'] = lamb.values.tolist()

    return config


def get_max_lamb(config):
    # Get layer names connected by bipartites
    positions_bipartites = [(config['bipartite'][bipartite]['source'],
                             config['bipartite'][bipartite]['target'])
                            for bipartite in config['bipartite']]

    # Create an empty dataframe with layer names as index and columns
    # that we'll fill where it's possible according to the bipartites
    max_lamb = pd.DataFrame(np.zeros((3, 3)),
                            index=config['multiplex'].keys(),
                            columns=config['multiplex'].keys())

    # fill each position corresponding to bipartite options
    for position in positions_bipartites:
        print(position[0], '<-->', position[1])
        max_lamb.loc[position] = 1
    # Future : Since we can inverse source and target layers)
    # could be conditionned by bipartite directionality
    max_lamb += max_lamb.transpose()
    # filling the diagonal to allow intra-layer exploration
    max_lamb += np.eye(3).astype(int)
    # normalise per rows
    max_lamb = max_lamb.div(max_lamb.sum(axis=1), axis=0)

    return max_lamb


def check_lamb(lamb, config):
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
    max_lamb = (get_max_lamb(config) > 0).astype(int)

    # Count number of locations where lamb has no non-zero value
    # and max_lamb is zero
    X = np.sum(np.sum((max_lamb - lamb) < 0))

    # If X > 0, lamb is not a valid probability matrix according to bipartites.
    if X > 0:
        return False
    # Else, lamb is a valid probability matrix according to bipartites.
    else:
        return True


def draw_lamb(df):
    """Draw the lamb matrix as a directed graph (networkx.DiGraph).
    Parameters
    ----------
    df: pandas.DataFrame
        The lamb matrix.
    e.g.:
    df = pd.DataFrame([[0.5, 0.5, 0],
                       [0.5, 0, 0.5],
                       [0, 0.5, 0.5]],
                       index = ['TF', 'Gene', 'Peak'],
                       columns = ['TF', 'Gene', 'Peak'])
    Returns
    -------
    None

    """
    fig, ax = plt.subplots()

    G = nx.from_pandas_adjacency(df-np.eye(len(df))*df,
                                 create_using=nx.DiGraph())
    pos = {list(G.nodes)[-i-1]: np.array((0, i))
           for i in range(-len(G.nodes), 0)}
    print(pos)

    nx.draw_networkx(G, with_labels=True, pos=pos,
                     node_size=1500, width=4, alpha=0.8, font_weight="bold",
                     arrows=True, connectionstyle='arc3, rad = 0.8')

    edge_labels = nx.get_edge_attributes(G, 'weight')

    my_draw_networkx_edge_labels(G,
                                 pos,
                                 edge_labels=edge_labels,
                                 font_color='k',
                                 font_size=12,
                                 label_pos=15,
                                 rad=0.8,
                                 rotate=False)
    plt.show()


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
        rad=0):

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
        Vertical align. {'center','top','bottom','baseline','center_baseline'}

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
        rotation_matrix = np.array([(0, 1), (-1, 0)])
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
