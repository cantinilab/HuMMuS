import os
import pytest
import shutil
import pickle
import pandas
import unittest
from pytest_dryrun import dryrun

from hummuspy.utils import logged
import hummuspy as hummus
import hummuspy.config

# activate dryrun passing the flag --dryrun in command line

class TestHumusPy(unittest.TestCase):

    def setUp(self):
        # path for test data
        test_dir = os.path.dirname(os.path.realpath(__file__))
        self.data_dir = os.path.join(test_dir, 'data', 'processed')
        self.out_dir = os.path.join(test_dir, 'data', 'chen_multilayer')
        self.tmp = os.path.join(test_dir, 'data', 'tmp_results')
        self.dryrun_dir = os.path.join(test_dir, 'data', 'dryrun_files')

    def tearDown(self):
        bipartite_file = os.path.join(self.out_dir, "bipartite", "rna_drugs.tsv")
        if os.path.exists(bipartite_file):
            os.remove(bipartite_file)
        
        multiplex_file = os.path.join(self.out_dir, "multiplex", "drug_drug_network.tsv")
        if os.path.exists(multiplex_file):
            os.remove(multiplex_file)
        
        config_extended_file = os.path.join(self.out_dir, "config", "config_extended.yml")
        if os.path.exists(config_extended_file):
            os.remove(config_extended_file)
        
        for file in os.listdir(self.tmp):
            path = os.path.join(self.tmp, file)
            os.remove(path)
            
    @dryrun
    def test_generalConfig_dry(self):
        config = {
            'multiplex': {
                'RNA': {
                    'graph_type': ['00'], 
                    'layers': ['multiplex/RNA/GENIE3.tsv']
                },
                'peaks': {
                    'graph_type': ['00'], 
                    'layers': ['multiplex/TF/PPI.tsvmultiplex/peaks/Cicero.tsv']
                }
            },
            'bipartite': {
                'bipartite/atac_rna.tsv': {
                    'graph_type': '00',
                    'source': 'RNA',
                    'target': 'peaks'
                }
            },
           'seeds': "seeds/seeds.txt",
           'self_loops': 0,
           'r': 0.7
        }
        keys = ['multiplex', 'bipartite', 'seeds', 'self_loops', 'r']
               
        assert len(config.keys()) == 5
        assert list(config.keys()) == keys
        
    def test_generalConfig(self):
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
        keys = ['multiplex', 'bipartite', 'seeds', 'self_loops', 'r']
        
        assert len(config.keys()) == 5
        assert list(config.keys()) == keys
            
    @dryrun
    def test_openConfig_dry(self):
        cpath = os.path.join( self.dryrun_dir, 'processed_config.pkl' )
        config = pickle.load( open(cpath, 'rb') )
        keys = ['bipartite', 'eta', 'lamb', 'multiplex', 'restart_proba', 'seeds', 'self_loops']
               
        assert len(config.keys()) == 7
        assert list(config.keys()) == keys
        assert type(config['lamb']) == pandas.DataFrame
        assert type(config['eta']) == pandas.Series
        assert len(config['eta']) == 3
        assert len(config['lamb']) == 3
        assert len(config['lamb'].iloc[0,]) == 3
        
    def test_openConfig(self):
        cpath = os.path.join( self.out_dir, "config", "grn_config.yml")
        config = hummus.config.open_config(cpath)
        keys = ['bipartite', 'eta', 'lamb', 'multiplex', 'restart_proba', 'seeds', 'self_loops']
               
        assert len(config.keys()) == 7
        assert list(config.keys()) == keys
        assert type(config['lamb']) == pandas.DataFrame
        assert type(config['eta']) == pandas.Series
        assert len(config['eta']) == 3
        assert len(config['lamb']) == 3
        assert len(config['lamb'].iloc[0,]) == 3
            
    @dryrun
    def test_getLambPlot_dry(self):
        plpatha = os.path.join( self.dryrun_dir, "grn_plot.png" )
        plpathb = os.path.join( self.dryrun_dir, "enhancer_plot.png" )
        plpathc = os.path.join( self.dryrun_dir, "br_plot.png" )
        plpathd = os.path.join( self.dryrun_dir, "tg_plot.png" )

        configExt_name = os.path.join( self.dryrun_dir, "extended_config.yaml") 

        assert os.path.isfile(configExt_name) == True
        assert os.path.isfile(plpatha) == True
        assert os.path.isfile(plpathb) == True
        assert os.path.isfile(plpathc) == True
        assert os.path.isfile(plpathd) == True
        
    def test_getLambPlot(self):
        cpath = os.path.join( self.out_dir, "config", "grn_config.yml")
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
        config['eta'] = pandas.Series([1, 0, 0, 0], index=['TF', 'peaks', 'RNA', 'drugs'])
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

        plpatha = os.path.join( self.tmp, "grn_plot.png" )
        config['lamb'] = hummus.config.get_grn_lamb(config, draw=False, save_plot=True, name_plot=plpatha )
        plpathb = os.path.join( self.tmp, "enhancer_plot.png" )
        config['lamb'] = hummus.config.get_enhancers_lamb(config, draw=True, save_plot=True, name_plot=plpathb )
        plpathc = os.path.join( self.tmp, "br_plot.png" )
        config['lamb'] = hummus.config.get_binding_regions_lamb(config, draw=True, save_plot=True, name_plot=plpathc )
        plpathd = os.path.join( self.tmp, "tg_plot.png" )
        config['lamb'] = hummus.config.get_target_genes_lamb(config, draw=True, save_plot=True, name_plot=plpathd )

        configExt_name = os.path.join( self.tmp, "extended_config.yaml") # the beautiful name you chose for your config
        hummus.config.save_config(config, configExt_name )

        assert os.path.isfile(configExt_name) == True
        assert os.path.isfile(plpatha) == True
        assert os.path.isfile(plpathb) == True
        assert os.path.isfile(plpathc) == True
        assert os.path.isfile(plpathd) == True
            
    @dryrun
    def test_defineGrn_dry(self):
        out_path = os.path.join( self.dryrun_dir, "ranked_grn_out.tsv" )
        df = pandas.read_csv( out_path, sep='\t' )
        cnf = {}
        result = (df, cnf)
        columns = ['index', 'layer', 'path_layer', 'score', 'gene', 'tf']

        assert os.path.isfile(out_path) == True
        assert type(result) == tuple
        assert len(df.columns) == 6
        assert list(df.columns) == columns
        assert len(df) == 62000
        
    def test_defineGrn(self):
        cpath = os.path.join( self.out_dir, "config", "grn_config.yml")
        config = hummus.config.open_config(cpath)
        multilayer_folder = "./data/chen_multilayer"

        save_flag = True
        out_path = os.path.join( self.tmp, "ranked_grn_out.tsv" )
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
        
        columns = ['index', 'layer', 'path_layer', 'score', 'gene', 'tf']
        df, config = result
        
        assert os.path.isfile(out_path) == True
        assert type(result) == tuple
        assert len(df.columns) == 6
        assert list(df.columns) == columns
        assert len(df) == 62000
        
