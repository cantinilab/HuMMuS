import os
import pandas as pd

def _process_dgidb(in_dir, out_dir):
    path = os.path.join( in_dir, 'dgidb_drug_gene_interactions.tsv' )
    df = pd.read_csv( path, sep='\t')
    resdf = df[ ['gene_name', 'drug_claim_name', 'interaction_score'] ]
    opath = os.path.join(out_dir, 'rna_drugs.tsv')
    resdf.to_csv(opath, sep='\t', index=None, header=None)
    
def _process_ddinterdb(in_dir, out_dir):
    levels = { 'Moderate': 0.5, 'Minor': 0.1, 'Major': 1, 'Unknown': 0 }
    path = os.path.join( in_dir, 'ddinter_drug_drug_code_L.csv' )
    df = pd.read_csv( path, sep=',')
    resdf = df[ ['Drug_A', 'Drug_B', 'Level'] ]
    print(resdf.Level.unique())
    resdf['Level'] = resdf['Level'].apply( lambda x: levels[x] )
    opath = os.path.join(out_dir, 'drug_drug_network.tsv')
    resdf.to_csv(opath, sep='\t', index=None, header=None)

def process_raw_files():
    data_dir = os.path.dirname(os.path.realpath(__file__))
    in_dir = os.path.join(data_dir, 'raw')
    out_dir = os.path.join(data_dir, 'processed')
    
    _process_dgidb(in_dir, out_dir)
    _process_ddinterdb(in_dir, out_dir)


if (__name__ == "__main__" ):
    process_raw_files()
