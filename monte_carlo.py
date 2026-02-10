import pandas as pd
import numpy as np
import json
import networkx as nx
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from Filtro.standardizer import Standardizer
from os import makedirs

from utils import mol_normalizer, mol_similarity

def check_cofactors(prod, cofactors):
    avoid = []
    for p in prod:
        for c in cofactors:
            if mol_similarity(p, c, mode='smile') == 1:
                avoid.append(p)
                break
    new_prod = []
    for p in prod:
        if p not in avoid:
            new_prod.append(p)
    return new_prod
                
def generator(smile, r, filename, ec_num, leg_id):
    a=Standardizer()
    m=Chem.MolFromSmiles(smile)
    Chem.SanitizeMol(m)
    mol1 = a.sequence_rr_legacy(m)
    mol2 = Chem.AddHs(m)
    mol = [mol1, mol2]
    
    prod = []
    res = []
    for mm in mol:
        for i in range(len(r)):
            try:
                pr = r[i].RunReactants((mm,))
                if len(pr) > 0:
                    m = pr[0][0]
                    p = mol_normalizer(Chem.MolToSmiles(m), mode='smile')
                    for j in p.split('.'):
                        m = Chem.MolFromSmiles(j)
                        if Descriptors.MolWt(m) <= 1000:
                            prod.append(j)
                            res.append(j+'$'+smile+'$'+leg_id[i]+'$'+ec_num[i])
            except:
                continue
    res = list(dict.fromkeys(res))
    res = pd.DataFrame(res, columns=['Resultado'])
    res.to_csv(filename, mode='a')
    return prod

def prod_clean(prod):
    for i in range(len(prod)):
        a = prod[i].split('.')
        prod[i] = a[0]
        for i in a[1:]:
            prod.append(i)
    prod = list(dict.fromkeys(prod))
    return prod

def next_step(p_inicial, best_p, smile, r, p_objective):
    for i in p_inicial:
        smile.add(i)
    p_final=generator(p_inicial[best_p],r)
    p_final = prod_clean(p_final)
    p = []
    for i in p_final:
        if i not in smile:
            p.append(i)
    p_final = p
    matrix = np.zeros((len(p_final),1))
    for i in range(len(p_final)):
        matrix[i] = mol_similarity(p_final[i], p_objective, mode='smile')
    return p_final, matrix

def create_nodes(G, prod, p_objective, father):
    G.add_node(prod, score=mol_similarity(prod, p_objective, mode='smile'), simulation=1)
    G.add_edge(father, prod)

def save_graph(G, filename):
    data = nx.node_link_data(G)
    with open(filename, 'w') as f:
        f.write(json.dumps(data))
    
def select_node(G, inicial, total_simulation):
    childs = list(G.successors(inicial))
    a = nx.get_node_attributes(G, 'score')
    b = nx.get_node_attributes(G, 'simulation')
    while len(childs) != 0:
        w_childs = []
        for i in childs:
            w_childs.append(a[i]+(np.sqrt(2)*np.sqrt(np.log(b[i])/total_simulation)))
        w_childs = np.array(w_childs)
        w_max_idx = w_childs.argmax()
        inicial = childs[w_max_idx]
        childs = list(G.successors(inicial))
        b[inicial] = 1
        nx.set_node_attributes(G, b, 'simulation')
    return inicial

def backpropagation(G, childs):
    father = list(G.predecessors(childs[0]))
    while len(father) != 0:
        childs = list(G.successors(father[0]))
        a = nx.get_node_attributes(G, 'score')
        w_childs = []
        for i in childs:
            w_childs.append(a[i])
        w_max = max(w_childs)
        a[father[0]] = w_max
        nx.set_node_attributes(G, a, 'score')
        father = list(G.predecessors(father[0]))    
    
def no_childs(G, father):
    a = nx.get_node_attributes(G, 'score')
    a[father] = 0.0
    nx.set_node_attributes(G, a, 'score')
    father = list(G.predecessors(father))
    while len(father) != 0:
        childs = list(G.successors(father[0]))
        a = nx.get_node_attributes(G, 'score')
        w_childs = []
        for i in childs:
            w_childs.append(a[i])
        w_max = max(w_childs)
        a[father[0]] = w_max
        nx.set_node_attributes(G, a, 'score')
        father = list(G.predecessors(father[0]))

def count_simulation(G):
    a = nx.get_node_attributes(G, 'simulation')
    for i in a:
        a[i] += 1
    nx.set_node_attributes(G, a, 'simulation')

inicial = pd.read_csv('data/Detectable.csv',usecols=['SMILES'])
inicial = inicial['SMILES'].tolist()

producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
producible = list(dict.fromkeys(producible['SMILES']))

cofactors = pd.read_csv('data/cofactors.csv', usecols=['Smiles'])
cofactors = cofactors['Smiles'].tolist()

DETECTABLE_ID = 0
PRODUCIBLE_ID = 0

start = inicial[DETECTABLE_ID]
smile = set([start])
end = producible[PRODUCIBLE_ID]

makedirs('generator_results/D'+str(DETECTABLE_ID), exist_ok=True)

filename = 'generator_results/D'+str(DETECTABLE_ID)+'/Results_D'+str(DETECTABLE_ID)+'P'+str(PRODUCIBLE_ID)+'.csv'
graph_filename = 'generator_results/D'+str(DETECTABLE_ID)+'/Graph_D'+str(DETECTABLE_ID)+'P'+str(PRODUCIBLE_ID)+'.json'

if mol_similarity(start, end, mode='smile') == 1:
    print('Detectable and producible are the same')
    exit()

rule_file = 'data/reactions_d16.csv'
df=pd.read_csv(rule_file ,usecols=['Rule','EC Number','MNX ID', 'Usage'])
ec_num=df["EC Number"].to_list()
leg_id=df["MNX ID"].to_list()
rule=df["Rule"].tolist()
r=[]
ec_num = []
leg_id = []
for i in range(len(rule)):
    if df['Usage'][i] != 'forward':
        r.append(Chem.rdChemReactions.ReactionFromSmarts(rule[i]))
        ec_num.append(df['EC Number'][i])
        leg_id.append(df['MNX ID'][i])
del (rule, df)

try:
    with open(graph_filename) as f:
        data = json.load(f)
    G = nx.node_link_graph(data)
    a = nx.get_node_attributes(G, 'simulation')
    val = []
    for i in a:
        val.append(a[i])
    total_simulation = max(val)
    del (val, a)
    for i in G.nodes():
        smile.add(i)
except:
    G = nx.DiGraph()
    G.add_node(start, score=mol_similarity(start, end, mode='smile'), simulation=1)
    total_simulation = 0

runs = 0
RDLogger.DisableLog('rdApp.*')
while nx.get_node_attributes(G, 'score')[start] != 1:
    runs += 1
    total_simulation += 1
    father = select_node(G, start, total_simulation)
    p = generator(father, r, filename, ec_num, leg_id)
    p = prod_clean(p)
    p = check_cofactors(p, cofactors)
    prod = []
    for i in p:
        if i not in smile:
            prod.append(i)
        smile.add(i)
    for i in prod:
        create_nodes(G, i, end, father)
    count_simulation(G)
    if len(prod) != 0:
        backpropagation(G, prod)
        print('Iteration done')
    else:
        no_childs(G, father)
        print('Iteration done, no childs')
    if nx.get_node_attributes(G, 'score')[start] == 0:
        print('Pathway not existant')
        break
    if runs == 5:
        save_graph(G, graph_filename)
        runs = 0

print('Script terminated')
save_graph(G, graph_filename)
