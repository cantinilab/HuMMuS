## Example of new bipartite and multiplex data processing

#### Gene-Drug interaction (bipartite - bipartite/rna_drugs.tsv)
The dataset in raw/dgidb\_drug\_gene\_interactions.tsv was downloaded from https://dgidb.org/downloads (accessed in Feb, 2026) release Dec, 2024. The pre-processing method extracts the column with gene names (gene\_name), drug names (drug\_name) and interaction score (interaction\_score - edge weight) and exports a tab-separated table file with these three columns without index and header.

#### Drug-Drug interaction (multiplex for drugs - multiplex/drug_drug_network.tsv)
The dataset in raw/ddinter\_drug\_drug\_code\_L.csv was downloaded from https://ddinter.scbdd.com/download/ (accessed in Feb, 2026), it only contains the interactions among drugs belonging to the ATC L category (antineoplastic and immunomodulating agents ). The pre-processing method extracts the columns with the drug partners (Drug\_A, Drug\_B) and uses the Level column to add the edge weight. The level is given as a categorical variable, whose values may be Major, Moderate, Minor or Unknown that were respectively mapped to 1, 0.5, 0.2 and 0. The pre-processing script also exports a tab-separated table file with these three columns without index and header.

