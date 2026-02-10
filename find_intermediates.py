import json
import pandas as pd
import networkx as nx

from os import listdir
from utils import mol_normalizer, rule_id_converter, smile_in_sink

pair_name = 'D0P27'
path_formated = pd.read_csv('MC/'+pair_name+'/MC_Pathways.csv')
cmpd_formated = pd.read_csv('MC/'+pair_name+'/MC_Compounds.tsv', sep='\t')
res_formated = pd.read_csv('MC/'+pair_name+'/MC_Scope.csv')

ecoli = pd.read_csv('data/ECOLI.csv')
ecoli = list(dict.fromkeys(ecoli['SMILES']))
cofactor = pd.read_csv('data/cofactors.csv')
cofactor = list(dict.fromkeys(cofactor['Smiles']))
producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
producible = list(dict.fromkeys(producible['SMILES']))

sink = ecoli + cofactor + producible
sink = list(dict.fromkeys(sink))

right = path_formated['Right'].tolist()
right = list(dict.fromkeys(right))
for i in range(len(right)):
    right[i] = right[i][2:]

compounds = cmpd_formated['Compound ID'].tolist()
structure = cmpd_formated['Structure'].tolist()

intermedia = []
for i in range(len(compounds)):
    try:
        if compounds[i] not in right and not smile_in_sink(i, sink):
            intermedia.append(structure[i])
    except:
        continue