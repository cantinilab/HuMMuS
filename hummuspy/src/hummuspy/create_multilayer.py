import networkx
import numpy
import os
import pandas
import sys

from multixrank.logger_setup import logger
import itertools
import sys

import numpy
import pandas

from multixrank import logger_setup
from multixrank import logger_setup
from multixrank.MultiplexAll import MultiplexAll
from multixrank import logger_setup, constants
from multixrank.BipartiteAll import BipartiteAll
from multixrank.MultiplexAll import MultiplexAll
from multixrank.Multiplex import Multiplex
from multixrank.ParameterEta import ParameterEta
from multixrank.ParameterLambda import ParameterLambda
from multixrank.Parameters import Parameters
from multixrank.logger_setup import logger
from multixrank.TransitionMatrix import TransitionMatrix

from typing import Union

class Bipartite:

    """Bipartite layer"""

    def __init__(self, key, abspath, graph_type, self_loops, on_disk=True, edge_list_df=None):

        """

        Args:
            abspath: str
            existing absolute path

            graph_type: str
            takes values 00=(unweighted, undirected), 01=(unweighted, directed),
            10=(weighted, undirected), 11=(weighted, directed)
        """

        self.key = key
        self.abspath = abspath
        self.graph_type = graph_type
        self.self_loops = self_loops
        self.edge_list_df = edge_list_df
        self.on_disk = on_disk

        if self.on_disk == True:
            if not os.path.isfile(abspath):  # error if path not exist
                logger.error("This path does not exist: {}".format(abspath))
                sys.exit(1)

        if not (graph_type in ['00', '10', '01', '11']):
            logger.error('MultiplexLayer multigraph type must take one of these values: 00, 10, 01, 11. '
                         'Current value: {}'.format(graph_type))
            sys.exit(1)
        self.graph_type = graph_type
        self.self_loops = self_loops

        self._networkx = None

    @property
    def networkx(self) -> networkx.Graph:
        """Converts layer to multigraph networkx object"""

        if self._networkx is None:

            names = ['col2', 'col1']  # layer file column labels changed
            dtype = str
            edge_attr = ['network_key']
            usecols = [0, 1]  # two cols like in unweighted
            if self.graph_type[1] == '1':  # unweighted/weighted layer (0/1)
                names = ['col2', 'col1', 'weight'] # changed
                dtype = {'col1': str, 'col2': str, 'weight': numpy.float64}
                edge_attr = ['network_key', 'weight']
                usecols = [0, 1, 2]  # two cols like in unweighted

            networkx_graph_obj = networkx.Graph()  # layer file column labels
            if self.graph_type[0] == '1':  # undirected/directed layer (0/1)
                networkx_graph_obj = networkx.DiGraph()

            if self.on_disk == True:
                multiplex_layer_edge_list_df = pandas.read_csv(self.abspath, sep="\t", header=None, names=names, dtype=dtype, usecols=usecols)
            else:
                multiplex_layer_edge_list_df = self.edge_list_df

            # remove df lines with self-loops, ie source==target
            if not self.self_loops:
                multiplex_layer_edge_list_df = multiplex_layer_edge_list_df.loc[
                    ~(multiplex_layer_edge_list_df.col1 == multiplex_layer_edge_list_df.col2)]
            multiplex_layer_edge_list_df['network_key'] = self.key

            self._networkx = networkx.from_pandas_edgelist(
                df=multiplex_layer_edge_list_df, source='col2', target='col1',
                edge_attr=edge_attr, create_using=networkx_graph_obj) # changed

            self._networkx.remove_edges_from(networkx.selfloop_edges(self._networkx))

            # networkx has no edges
            # TODO replace edges with nodes
            if len(self._networkx.edges()) == 0:
                logger.error(
                    'The following bipartite graph does not return any edge: {}'.format(
                        self.key))
                sys.exit(1)

        return self._networkx


class MultiplexLayer:

    """Multiplex layer"""

    def __init__(self, key, abspath, graph_type, multiplex, tau, self_loops, on_disk=True, edge_list_df=None):

        """

        Args:
            abspath: str
            existing absolute path

            graph_type: str
            takes values 00=(unweighted, undirected), 01=(unweighted, directed),
            10=(weighted, undirected), 11=(weighted, directed)

            multiplex: str
            Parent multiplex key
        """

        self.key = key
        self.abspath = abspath
        self.multiplex = multiplex
        self.tau = tau
        self.graph_type = graph_type
        self.self_loops = self_loops
        self.edge_list_df = edge_list_df
        self.on_disk = on_disk

        if self.on_disk == True:
            if not os.path.isfile(abspath):  # error if path not exist
                logger.error("This path does not exist: {}".format(abspath))
                sys.exit(1)

        if not (graph_type in ['00', '10', '01', '11']):
            logger.error('MultiplexLayer multigraph type must take one of these values: 00, 10, 01, 11. '
                         'Current value: {}'.format(graph_type))
            sys.exit(1)

        self._networkx = None

    @property
    def networkx(self) -> networkx.Graph:
        """Converts layer to multigraph networkx object"""

        if self._networkx is None:

            # if multigraph is None:
            name_lst = ['source', 'target']  # layer file column labels
            dtype_dic = dict(zip(name_lst, [str]*len(name_lst)))

            # unweighted vs weighted ######################################
            edge_attr = ['network_key']
            usecols = [0, 1]  # two cols like in unweighted
            if self.graph_type[1] == '1':  # unweighted/weighted layer (0/1)
                name_lst = ['source', 'target', 'weight']
                dtype_dic['weight'] = numpy.float64
                edge_attr = ['network_key', 'weight']
                usecols = [0, 1, 2]  # three cols in weighted

            # undirected vs directed ######################################
            networkx_graph_obj = networkx.Graph()  # layer file column labels
            if self.graph_type[0] == '1':  # undirected/directed layer (0/1)
                networkx_graph_obj = networkx.DiGraph()

            if self.on_disk == True:
                multiplex_layer_edge_list_df = pandas.read_csv(
                    self.abspath, sep="\t", header=None, names=name_lst, dtype=dtype_dic, keep_default_na=False, usecols=usecols)
            else:
                multiplex_layer_edge_list_df = self.edge_list_df
            # remove df lines with self-loops, ie source==target if bipartite_notes=true
            if not self.self_loops:
                multiplex_layer_edge_list_df = multiplex_layer_edge_list_df.loc[
                    ~(multiplex_layer_edge_list_df.source == multiplex_layer_edge_list_df.target)]
            multiplex_layer_edge_list_df['network_key'] = self.key
            self._networkx = networkx.from_pandas_edgelist(
                df=multiplex_layer_edge_list_df, source='source',
                target='target', create_using=networkx_graph_obj, edge_attr=edge_attr)

            self._networkx.remove_edges_from(networkx.selfloop_edges(self._networkx))

            # networkx has no edges
            if len(self._networkx.nodes()) == 0:
                logger.error(
                    'The following multiplex layer does not return any edge: {}'.format(self.key))
                sys.exit(1)

        return self._networkx


class Seed:

    """Handle seeds"""

    def __init__(self, path: str, multiplexall: MultiplexAll, seeds: list=list()):
        """Handles seeds

        Args:
            path (str): seed path
            node (list): List of nodes in multiplexes

        Returns:
            none

        """

        # self._seed_df = pandas.read_csv(path, sep="\t", header=None, names=["seeds"], dtype=str)
        self._path = path
        self._multiplexall_obj = multiplexall

        self._seed_list = seeds
        self._multiplex_seed_list2d = list()

    @property
    def seed_list(self) -> list:
        """Get list with all nodes"""


        ###################################################################
        #
        # Check seed list
        #
        ###################################################################

        if len(self._seed_list) <= 0:  # no seeds
            logger_setup.logger.error(
                "These seed nodes were not found in the network: {}".format(
                    self._seed_list))
            sys.exit(1)

        seeds_not_in_multiplex_set = set(self._seed_list) - set(self._multiplexall_obj.nodes)
        if len(seeds_not_in_multiplex_set) > 0:  # seeds not in multiplex
            logger_setup.logger.error(
                "These seed nodes were not found in the network: {}".format(
                    seeds_not_in_multiplex_set))
            sys.exit(1)

        return self._seed_list

    @property
    def multiplex_seed_list2d(self) -> list:
        """Returns valid seed label list"""

        if not self._multiplex_seed_list2d:

            for multiplex in self._multiplexall_obj.multiplex_tuple:

                self._multiplex_seed_list2d.append(list(set(self._seed_list) & set(multiplex.nodes)))

        return self._multiplex_seed_list2d

    # 1.3.1 :
    def get_seed_scores(self, transition) -> (numpy.array, pandas.DataFrame):
        """

        Function that determine the initial probability distribution thanks to seeds and normalization.
        For the seed_rank() and homogeneous_seed_rank() functions.

        Returns :
                seed_score (pandas.core.frame.DataFrame) : A Dataframe with value of probaility for the seeds in each layer
                The value of probability is zero outhere from seeds.

        """

        multiplexall_layer_count_list = [len(m.layer_tuple) for m in self._multiplexall_obj.multiplex_tuple]
        multiplexall_node_count_list = [len(m.nodes) for m in self._multiplexall_obj.multiplex_tuple]
        multiplexall_node_list2d = [multiplexone_obj.nodes for multiplexone_obj in self._multiplexall_obj.multiplex_tuple]

        multiplexall_layer_key_list = []
        multiplexall_layer_key_list2d = []
        for multiplex in self._multiplexall_obj.multiplex_tuple:
            multiplexall_layer_key_list2d.append([layer_obj.key for layer_obj in multiplex.layer_tuple])
            for layer in multiplex.layer_tuple:
                multiplexall_layer_key_list.append(layer.key)

        multiplexall_layer_key_lst = [layer.key for multiplex in self._multiplexall_obj.multiplex_tuple for layer in multiplex.layer_tuple]
        seed_score_df = pandas.DataFrame(0, index=self._seed_list, columns=multiplexall_layer_key_lst)
        prox_vector = numpy.zeros((numpy.shape(transition)[0], 1))
        for seed_label in seed_score_df.index:  # loop through seeds
            for multiplex_idx, multiplex in enumerate(self._multiplexall_obj.multiplex_tuple):
                if seed_label in multiplex.nodes:
                    for layer_idx, layer in enumerate(multiplex.layer_tuple):
                        seed_score_df.loc[seed_label, layer.key] = multiplex.eta * layer.tau / len(self.multiplex_seed_list2d[multiplex_idx])
                        start = sum(numpy.array(multiplexall_node_count_list[:multiplex_idx]) * numpy.array(multiplexall_layer_count_list[:multiplex_idx])) + (multiplexall_node_count_list[multiplex_idx] * layer_idx)
                        pos = multiplexall_node_list2d[multiplex_idx].index(seed_label) + start
                        prox_vector[pos] = seed_score_df.loc[seed_label, layer.key]
        return prox_vector, seed_score_df


class Multixrank:
    """Main class to run the random walk with restart in universal multiplex networks"""

    def __init__(
        self,
        multiplex,
        bipartite,
        eta,
        lamb,
        seeds,
        self_loops,
        restart_proba,
        pr=None):
        """
        Constructs an object for the random walk with restart.

        Args:
            config (str): Path to the configuration file in YML format. Paths will be used relative to the wdir path variable below
            wdir (str): Path to the working directory that will be as starting point to the paths in the config file.
        """

        #######################################################################
        #
        # Read ConfigPath
        #
        #######################################################################

        config_parser_obj = CreateMultilayer(
            multiplex=multiplex,
            bipartite=bipartite,
            eta=eta,
            lamb=lamb,
            seeds=seeds,
            self_loops=self_loops,
            restart_proba=restart_proba)
            
        config_parser_obj.parse()
        
        self.pr = pr

        #######################################################################
        #
        # paramater object from config_parser and properties
        #
        #######################################################################
        parameter_obj = config_parser_obj.parameter_obj
        self.r = parameter_obj.r
        self.lamb = parameter_obj.lamb
        logger.debug("Parameter 'lambda' is equal to: {}".format(self.lamb))

        #######################################################################
        #
        # multiplexall object from config_parser and properties
        #
        #######################################################################

        multiplexall_obj = config_parser_obj.multiplexall_obj
        self.multiplexall_obj = multiplexall_obj

        #######################################################################
        #
        # bipartite object from config_parser and properties
        #
        #######################################################################

        self.bipartiteall_obj = config_parser_obj.bipartitelist_obj

        self.multiplex_layer_count_list = [len(multiplexone_obj.layer_tuple) for multiplexone_obj in multiplexall_obj.multiplex_tuple]

        self.multiplexall_node_list2d = [multiplexone_obj.nodes for multiplexone_obj in multiplexall_obj.multiplex_tuple]

        # self. N nb of nodes in each multiplex
        # self.N = list()
        self.multiplexall_node_count_list = [len(x) for x in self.multiplexall_node_list2d]

        #######################################################################
        #
        # seed object from config_parser and properties
        #
        #######################################################################
        if type(self.pr) == type(None):
            self.seed_obj = config_parser_obj.seed_obj
        else:
            N = copy.deepcopy(self.multiplexall_node_count_list)
            N.insert(0,0)
            L = self.multiplex_layer_count_list
            temp = [numpy.repeat((self.pr[numpy.sum(N[:i]):numpy.sum(N[:i+1])]/L[i-1]),L[i-1]) for i in range(1,len(L)+1)]
            self.pr = numpy.concatenate(temp)

    # 1.3.3 :
    def __random_walk_restart(self, prox_vector, transition_matrixcoo, r):
        """

        Function that realize the RWR and give back the steady probability distribution for
        each multiplex in a dataframe.

        self.results (list) : A list of ndarray. Each ndarray correspond to the probability distribution of the
            nodes of the multiplex.

        """
        rwr_result_lst = list()
        threshold = 1e-10
        residue = 1
        itera = 1
        prox_vector_norm = prox_vector / (sum(prox_vector))
        restart_vector = prox_vector_norm
        while residue >= threshold:
            old_prox_vector = prox_vector_norm
            prox_vector_norm = (1 - r) * (transition_matrixcoo.dot(prox_vector_norm)) + r * restart_vector
            residue = numpy.sqrt(sum((prox_vector_norm - old_prox_vector) ** 2))
            itera += 1
        for k in range(len(self.multiplex_layer_count_list)):
            start = sum(numpy.array(self.multiplexall_node_count_list[:k]) * numpy.array(self.multiplex_layer_count_list[:k]))
            end = start + self.multiplexall_node_count_list[k] * self.multiplex_layer_count_list[k]
            data = numpy.array(prox_vector_norm[start:end])
            rwr_result_lst.append(data)
        return rwr_result_lst

    ###########################################################################
    # 2 :Analysis func##############################tions
    ###########################################################################

    # 2.1 :
    ###########################################################################

    # 2.1.1 :
    def random_walk_rank(self) -> pandas.DataFrame:
        """
        Function that carries ous the full random walk with restart from a list of seeds.

        Returns :
                rwr_ranking_df (pandas.DataFrame) : A pandas Dataframe with columns: multiplex, node, layer, score
        """
        bipartite_matrix = self.bipartiteall_obj.bipartite_matrix
        transition_matrix_obj = TransitionMatrix(multiplex_all=self.multiplexall_obj, bipartite_matrix=bipartite_matrix, lamb=self.lamb)
        transition_matrixcoo = transition_matrix_obj.transition_matrixcoo

        # Get initial seed probability distribution
        if type(self.pr) == type(None):
            prox_vector, seed_score = self.seed_obj.get_seed_scores(transition=transition_matrixcoo)
        else:
            prox_vector = self.pr
        # Run RWR algorithm
        rwr_ranking_lst = self.__random_walk_restart(prox_vector, transition_matrixcoo, self.r)
        rwr_ranking_df = self.__random_walk_rank_lst_to_df(rwr_result_lst=rwr_ranking_lst)

        return rwr_ranking_df

    def write_ranking(self, random_walk_rank: pandas.DataFrame, path: str, top: int = None, aggregation: str = "gmean", degree: bool = False):
        """Writes the 'random walk results' to a subnetwork with the 'top' nodes as a SIF format (See Cytoscape documentation)

        Args:
            rwr_ranking_df (pandas.DataFrame) : A pandas Dataframe with columns: multiplex, node, layer, score, which is the output of the random_walk_rank function
            path (str): Path to the SIF file
            top (int): Top nodes based on the random walk score to be included in the SIF file
            aggregation (str): One of "nomean", "gmean", "hmean", "mean", or "sum"
        """

        if not (aggregation in ['nomean', 'gmean', 'hmean', 'mean', 'sum']):
            logger.error('Aggregation parameter must take one of these values: "nomean", "gmean", "hmean", "mean", or "sum". '
                         'Current value: {}'.format(aggregation))
            sys.exit(1)
        
        output_obj = Output(random_walk_rank, self.multiplexall_obj, top=top, top_type="layered", aggregation=aggregation)
        output_obj.to_tsv(outdir=path, degree=degree)

    def to_sif(self, random_walk_rank: pandas.DataFrame, path: str, top: int = None, top_type: str = 'layered', aggregation: str = 'gmean'):
        """Writes the 'random walk results' to a subnetwork with the 'top' nodes as a SIF format (See Cytoscape documentation)

        Args:
            rwr_ranking_df (pandas.DataFrame) : A pandas Dataframe with columns: multiplex, node, layer, score, which is the output of the random_walk_rank function
            path (str): Path to the TSV file with the random walk results
            top (int): Top nodes based on the random walk score to be included in the TSV file
            top_type (str): "per layer" (top nodes for each layer) or "all" (top nodes any layer)
            aggregation (str): One of "none", "geometric mean" or "sum"
        """

        if not (aggregation in ['nomean', 'gmean', 'hmean', 'mean', 'sum']):
            logger.error('Aggregation parameter must take one of these values: "nomean", "gmean", "hmean", "mean", or "sum". '
                         'Current value: {}'.format(aggregation))
            sys.exit(1)

        if not (top_type in ['layered', 'all']):
            logger.error('top_type parameter must take one of these values: "layered" or "all". '
                         'Current value: {}'.format(top_type))
            sys.exit(1)
        
        output_obj = Output(random_walk_rank, self.multiplexall_obj, top=top, top_type=top_type, aggregation=aggregation)
        pathlib.Path(os.path.dirname(path)).mkdir(exist_ok=True, parents=True)
        output_obj.to_sif(path=path, bipartiteall=self.bipartiteall_obj)

    def __random_walk_rank_lst_to_df(self, rwr_result_lst) -> pandas.DataFrame:
        rwrrestart_df = pandas.DataFrame({'multiplex': [], 'node': [], 'layer': [], 'score': []})
        for i, multiplex in enumerate(self.multiplexall_obj.multiplex_tuple):
            multiplex_label_lst = [multiplex.key] * len(multiplex.nodes) * len(
                multiplex.layer_tuple)
            nodes = [item for subl in [multiplex.nodes] * len(multiplex.layer_tuple) for item in
                     subl]
            layer_lst = [item for subl in
                         [[layer.key] * len(multiplex.nodes) for layer in multiplex.layer_tuple] for
                         item in subl]
            if type(self.pr) == type(None):
                score = list(rwr_result_lst[i].T[0])
            else:
                score = list(rwr_result_lst[i].T)
            rwrrestart_df = pandas.concat([rwrrestart_df, pandas.DataFrame(
                {'multiplex': multiplex_label_lst, 'node': nodes, 'layer': layer_lst, 'score': score})], axis=0)
        return rwrrestart_df

class CreateMultilayer:

    def __init__(self, multiplex, bipartite, eta, lamb, seeds, self_loops, restart_proba):

        """Takes a confMultiplexLayerig_path of the config yaml file"""

        
        
        # default is 00 for all interaction types
        self.multiplex = multiplex
        self.bipartite = bipartite
        self.eta = eta
        self.lamb = lamb
        self.self_loops = self_loops
        self.seeds = seeds
        self.restart_proba = restart_proba


    def parse(self):

        """Parses the config yaml file and give fields to specialized functions"""

        #######################################################################
        #
        # Parses multiplex layers
        # and create MultiplexLayer, MultiplexOne and MultiplexAll objects
        #
        #######################################################################

        self.multiplexall_obj = self.__parse_multiplex()
        self.multiplex = "Done"
        #######################################################################
        #
        # Parses bipartite
        #
        #######################################################################

        self.bipartitelist_obj = self.__parse_bipartite()
        self.bipartite = "Done"

        #######################################################################
        #
        # Parses Seeds
        #
        #######################################################################

        self.seed_obj = self.__parse_seed()
        self.seeds = "Done"
        #######################################################################
        #
        # Parameter 'delta' check
        #
        #######################################################################

#        if 'delta' in self.config_dic:  # user defined (changed)
 #           delta_lst = [float(Fraction(delta)) for delta in self.config_dic['delta']] # changed
  #          Parameters.check_eta(delta_lst , len(self.multiplexall_obj.multiplex_tuple)) # changed
            
            
            
        ###################################################################
        #
        # Parameter 'lamb' check
        #
        ###################################################################

        lamb = numpy.array(self.lamb)
        Parameters.check_lamb(lamb , len(self.multiplexall_obj.multiplex_tuple)) # changed
                

        #######################################################################
        #
        # Parameter 'eta': Parse and set
        #
        #######################################################################
        if self.eta != None:
            eta_lst = [float(Fraction(eta)) for eta in self.eta]
            Parameters.check_eta(eta_lst , len(self.multiplexall_obj.multiplex_tuple)) # changed
        else:  # default eta
            n = len(self.multiplexall_obj.multiplex_tuple)
            alpha = [len(x) for x in self.seed_obj.multiplex_seed_list2d]
            # Generate eta default
            eta_lst = ParameterEta(n, alpha).vect_X().tolist()
        logger.debug("Parameter 'eta' is equal to: {}".format(eta_lst))

        for i,multiplex_obj in enumerate(self.multiplexall_obj.multiplex_tuple):
            multiplex_obj.eta = eta_lst[i]

        #######################################################################
        #
        # Parses parameters
        #
        #######################################################################

        seed_count_list2d = [len(i) for i in self.seed_obj.multiplex_seed_list2d]
        self.parameter_obj = self.__parse_parameters(seed_count_list2d=seed_count_list2d)

    def __parse_bipartite(self):
        """
        Reads multiplex field and create MultiplexAll object
        """

        source_target_bipartite_dic = {}

        # loop over each multiplex
        for i, layer_key in enumerate(self.bipartite):

            graph_type = '00'
            if 'graph_type' in self.bipartite[layer_key]:
                graph_type = self.bipartite[layer_key]['graph_type']

            if 'source' in self.bipartite[layer_key]:
                bipartite_source = self.bipartite[layer_key]['source']
            else:
                logger_setup.logger.error("No 'source' field found for bipartite network")
                sys.exit(1)

            if 'target' in self.bipartite[layer_key]:
                bipartite_target = self.bipartite[layer_key]['target']
            else:
                logger_setup.logger.error("No 'target' field found for bipartite network")
                sys.exit(1)

            layer_obj = Bipartite(
                key=layer_key,
                abspath='',
                graph_type=graph_type,
                self_loops=self.self_loops,
                on_disk=False,
                edge_list_df = self.bipartite[layer_key]['edge_list_df'])
            source_target_bipartite_dic[(bipartite_source, bipartite_target)] = layer_obj

        return BipartiteAll(source_target_bipartite_dic, multiplexall=self.multiplexall_obj)

    def __parse_multiplex(self):
        """
        Reads multiplex field and create MultiplexAll object
        """


        multiplex_count = len(self.multiplex)
        # convert int to strings
        multiplex_obj_list = []

        # loop over each multiplex
        for multiplex_idx, multiplex_key in enumerate(self.multiplex):

            layer_obj_list = []
            multiplex_node_list = []
            layer_key_tuple = tuple(self.multiplex[multiplex_key]['names'])
            print(layer_key_tuple)
            ###################################################################
            #
            # tau
            #
            ###################################################################

            tau_lst = [float(Fraction(1/len(layer_key_tuple)))] * len(layer_key_tuple)
            if 'tau' in self.multiplex[multiplex_key]:
                tau_lst  = [float(Fraction(tau)) for tau in self.multiplex[multiplex_key]['tau']]
            Parameters.check_tau(tau_lst , len(layer_key_tuple))
            logger.debug("Multiplex '{}'. Parameter 'tau' is equal: {}".format(multiplex_key, tau_lst))


            ###################################################################
            #
            # layer multigraph types
            #
            ###################################################################

            graph_type_lst = ['00'] * len(layer_key_tuple)
            if 'graph_type' in self.multiplex[multiplex_key]:
                graph_type_lst = self.multiplex[multiplex_key]['graph_type']

            ###################################################################
            #
            # layers
            #
            ###################################################################

            # loop over layers
            for layer_idx, layer_key in enumerate(layer_key_tuple):
                print(layer_key)
                layer_obj = MultiplexLayer(key=layer_key,
                                           abspath='',
                                           graph_type=graph_type_lst[layer_idx],
                                           multiplex=multiplex_key,
                                           tau=tau_lst[layer_idx],
                                           self_loops=self.self_loops,
                                           on_disk=False,
                                           edge_list_df = self.multiplex[multiplex_key]['layers'][layer_idx])

                multiplex_node_list = sorted([*set(multiplex_node_list + [*layer_obj.networkx.nodes])])
                layer_obj_list.append(layer_obj)

            # Append missing nodes from other layers in one multiplex to each layer
            for layer_idx, layer_obj in enumerate(layer_obj_list):
                layer_missing_nodes = sorted(set(multiplex_node_list) - set(list(layer_obj.networkx.nodes)))
                layer_obj.networkx.add_nodes_from(layer_missing_nodes)
                layer_obj_list[layer_idx] = layer_obj  # update

            ###################################################################
            #
            # Parameter: delta
            # If only layer: delta=0
            # If more than one layer: delta=0.5
            #
            ###################################################################

            if 'delta' in self.multiplex[multiplex_key]:  # user-defined delta
                delta = float(self.multiplex[multiplex_key]['delta'])
            else:  # default delta
                if len(layer_obj_list) == 1:  # only one layer, delta=0
                    delta = 0
                else:  # more than one layer, delta=0.5
                    delta = 0.5
            logger.debug("Multiplex: {}. Parameter 'delta' is equal to: {}".format(multiplex_key, delta))

            ###################################################################
            #
            # multiplex object
            #
            ###################################################################

            multiplex_one_obj = Multiplex(multiplex_key, layer_tuple=tuple(layer_obj_list), delta=delta)
            multiplex_obj_list.append(multiplex_one_obj)

        return MultiplexAll(multiplex_tuple=tuple(multiplex_obj_list))

    def __parse_parameters(self, seed_count_list2d):

        r = constants.r
        if self.restart_proba != None:
            r = float(self.restart_proba)
        logger.debug("r is equal to: {}".format(r))
        
        lamb_arr = None
        lamb_2dlst = [[float(Fraction(j)) for j in i] for i in self.lamb]
        lamb_arr = numpy.array(lamb_2dlst)

        parameters_obj = Parameters(r=r, lamb=lamb_arr, multiplexall=self.multiplexall_obj, seed_count_list2d=seed_count_list2d)

        return parameters_obj

    def __parse_seed(self):

        seed_obj = Seed(path='', multiplexall=self.multiplexall_obj, seeds = self.seeds)

        return seed_obj
    

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
    ATAC_RNA_bipartite_graph_type: str = '00'):

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
        TF_layer_graph_type (Union[str, list], optional): Graph type for TF layer(s). 
            Defaults to '00'.
        ATAC_layer_graph_type (Union[str, list], optional): Graph type for ATAC layer(s).
            Defaults to '01'.
        RNA_layer_graph_type (Union[str, list], optional): Graph type for RNA layer(s).
            Defaults to '01'.
        TF_ATAC_bipartite_graph_type (str, optional): Graph type for TF-ATAC bipartite.
            Defaults to '00'.
        ATAC_RNA_bipartite_graph_type (str, optional): Graph type for ATAC-RNA bipartite.
            Defaults to '00'.

    Returns:
        dict: Formatted data as input for MultiXrank-HuMMuS class
    """
    if type(TF_layer) == pandas.DataFrame:
        TF_layer = [TF_layer]
    if type(ATAC_layer) == pandas.DataFrame:
        ATAC_layer = [ATAC_layer]
    if type(RNA_layer) == pandas.DataFrame:
        RNA_layer = [RNA_layer]
    
    if type(TF_layer_graph_type) == str:
        TF_layer_graph_type = [TF_layer_graph_type]*len(TF_layer)
    if type(ATAC_layer_graph_type) == str:
        ATAC_layer_graph_type = [ATAC_layer_graph_type]*len(ATAC_layer)
    if type(RNA_layer_graph_type) == str:
        RNA_layer_graph_type = [RNA_layer_graph_type]*len(RNA_layer)

    assert len(TF_layer) == len(TF_layer_graph_type), "TF_layer and TF_layer_graph_type must have the same length"
    assert len(ATAC_layer) == len(ATAC_layer_graph_type), "ATAC_layer and ATAC_layer_graph_type must have the same length"
    assert len(RNA_layer) == len(RNA_layer_graph_type), "RNA_layer and RNA_layer_graph_type must have the same length"

    keys_data = ['names', 'graph_type', 'layers']
    # Multiplex
    multiplex = {
        'TF': {key:[] for key in keys_data},
        'ATAC': {key:[] for key in keys_data},
        'RNA': {key:[] for key in keys_data}
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