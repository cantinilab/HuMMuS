import os
import pytest
import shutil
import unittest
from pytest_dryrun import dryrun

from hummuspy.utils import logged
import hummuspy.config

# activate dryrun passing the flag --dryrun in command line

class TestHumusPy(unittest.TestCase):

    def setUp(self):
        # path for test data
        test_dir = os.path.dirname(os.path.realpath(__file__))
        self.data_dir = os.path.join(test_dir, 'data', 'processed')
        self.out_dir = os.path.join(test_dir, 'data', 'chen_multilayer')
        print( self.out_dir )

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
            #shutil.rmtree(self.out_dir)
            
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
           'seed': "seeds/seeds.txt",
           'self_loops': 0,
           'r': 0.7
        }
               
        assert len(config.keys()) == 5
        
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
        assert len(config.keys()) == 5
    
    
        
