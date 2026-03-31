import json
import networkx as nx
import pandas as pd
import csv
import os
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from collections import Counter

from utils import mol_similarity

def scope_format(pair_name):
    d_idx = pair_name[pair_name.index('D')+1:pair_name.index('P')]
    p_idx = pair_name[pair_name.index('P')+1:]
    data = json.load(open('generator_results/Updated_graph/Graph_'+pair_name+'.json'))
    G = nx.node_link_graph(data)
    rule = pd.read_csv('data/reactions_d16.csv')
    rule = rule[rule['Usage'] != 'forward']
    rule.reset_index(inplace=True, drop=True)
    leg_id = rule['MNX ID'].tolist()
    score = rule['Score normalized'].tolist()
    rule_id = rule['Rule ID'].tolist()
    res = pd.read_csv('generator_results/D'+d_idx+'/Results_'+pair_name+'.csv')
    res=list(dict.fromkeys(res['Resultado']))
    try:
        res.remove('Resultado')
    except:
        res = res
    padre = []
    reac = []
    prod = []
    ec_num = []
    for i in res:
        a,b,c,d= i.split('$')
        if G.has_node(b):
            padre.append(b)
            reac.append(c)
            prod.append(a)
            ec_num.append(d)
    del (a,b,c,d)
    
    RDLogger.DisableLog('rdApp.*')
    padre_inchi = []
    prod_inchi = []
    for i in range(len(prod)):
        m = Chem.MolFromSmiles(prod[i])
        prod_inchi.append(Chem.MolToInchi(m))
        m = Chem.MolFromSmiles(padre[i])
        padre_inchi.append(Chem.MolToInchi(m))
    
    ite = [0]*len(padre)
    for i in range(len(prod)):
        ite[i] = int(nx.shortest_path_length(G, padre[i], padre[0])/2)
    
    res_score = []
    res_r_id = []
    res_r_smile = []
    for i in reac:
        idx = leg_id.index(i)
        res_score.append(score[idx])
        res_r_id.append(rule_id[idx])
        r = Chem.rdChemReactions.ReactionFromSmarts(rule['Rule'][idx])
        rule_smile = Chem.rdChemReactions.ReactionToSmiles(r)
        res_r_smile.append(rule_smile)
    
    trans_id = []
    ite_reac = []
    for i in range(len(padre)):
        ite_reac.append(padre[i]+'$'+res_r_id[i])
    ite_reac = list(dict.fromkeys(ite_reac))
    for i in range(len(ite_reac)):
        ite_reac[i] = ite_reac[i].split('$')
    ite_reac[0].append(0)
    for i in range(1,len(ite_reac)):
        if ite_reac[i][0] != ite_reac[i-1][0]:
            ite_reac[i].append(0)
        else:
            ite_reac[i].append(ite_reac[i-1][2]+1)
    
    for i in range(len(prod)):
        for j in range(len(ite_reac)):
            if ite_reac[j][0] == padre[i] and ite_reac[j][1] == res_r_id[i]:
                t = 'TRS_0_'+str(ite[i])+'_'+str(ite_reac[j][2])
                break
        trans_id.append(t)
    
    for i in range(len(ec_num)):
        ec_num[i] = '[' + ec_num[i] +']'
        res_r_id[i] = '[' + res_r_id[i] +']'
        
    res_formated = pd.DataFrame()
    res_formated['Initial source']=['[target]']*len(padre)
    res_formated['Transformation ID']=trans_id
    res_formated['Reaction SMILES']=res_r_smile
    res_formated['Substrate SMILES']=padre
    res_formated['Substrate InChI']=padre_inchi
    res_formated['Product SMILES']=prod
    res_formated['Product InChI']=prod_inchi
    res_formated['In Sink']=[0]*len(res_formated)
    res_formated['Sink name']=['[None]']*len(res_formated)
    res_formated['Diameter']=[16]*len(res_formated)
    res_formated['Rule ID']=res_r_id
    res_formated['EC number']=ec_num
    res_formated['Score']=res_score
    res_formated['Starting Source SMILES']=[padre[0]]*len(res_formated)
    res_formated['Iteration']=ite

    return res_formated
    
def path_cmpd_format(pair_name):
    d_idx = int(pair_name[pair_name.index('D')+1:pair_name.index('P')])
    p_idx = int(pair_name[pair_name.index('P')+1:])
    
    data = json.load(open('generator_results/Updated_graph/Graph_'+pair_name+'.json'))
    G = nx.node_link_graph(data)
    
    inicial = pd.read_csv('data/Detectable.csv',usecols=['SMILES'])
    inicial = inicial['SMILES'].tolist()
    
    producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
    producible = list(dict.fromkeys(producible['SMILES']))
    producible.remove(producible[-1])
    
    for i in G.nodes:
        if list(G.predecessors(i)) == [] and mol_similarity(producible[p_idx], i, mode='smile') == 1:
            source = i
    if G.has_node(producible[p_idx]):
        source = producible[p_idx]
    
    rule = pd.read_csv('data/reactions_d16.csv')
    rule = rule[rule['Usage'] != 'forward']
    rule.reset_index(inplace=True, drop=True)
    r_id = rule['Rule ID'].tolist()
    leg_id = rule['MNX ID'].tolist()
    
    all_paths = []
    for i in nx.shortest_simple_paths(G, source, inicial[d_idx]):
        all_paths.append(i)
    
    nodes = []
    names = []
    idx = 2
    for i in G.nodes:
        if i[0] != 'M':
            nodes.append(i)
            if i == inicial[d_idx]:
                names.append('TARGET_0000000001')
            else:
                names.append('CMPD_'+str(idx).zfill(10))
                idx += 1
            
    left = []
    right = []
    rule_id = []
    path_id = []
    for i in range(len(all_paths)):
        for j in range(len(all_paths[i])):
            if all_paths[i][j][0] == 'M':
                rule_idx = leg_id.index(all_paths[i][j])
                rule_id.append(r_id[rule_idx])
                path_id.append(i)
                pred = list(G.predecessors(all_paths[i][j]))
                l_name = ''
                for k in pred:
                    l_idx = nodes.index(k)
                    l_name +='1.'+ names[l_idx] + ':'
                if l_name[-1] == ':':
                    left.append(l_name[0:-1])
                else:
                    left.append(l_name)
                r_idx = nodes.index(all_paths[i][j+1])
                right.append('1.'+names[r_idx])
    
    all_lr = []
    for i in range(len(left)):
        all_lr.append(left[i]+'/'+right[i])
    all_lr = list(dict.fromkeys(all_lr))
    
    unique_id = []
    res_formated = pd.read_csv('MC/'+pair_name+'/MC_Scope.csv', usecols=['Transformation ID', 'Substrate SMILES','Rule ID'])
    for i in range(len(path_id)):
        idx = names.index(right[i][2:])
        for j in range(len(res_formated)):
            if '['+rule_id[i]+']' == res_formated['Rule ID'][j] and nodes[idx] == res_formated['Substrate SMILES'][j]: 
                trs = res_formated['Transformation ID'][j] +'_0'
                break
        unique_id.append(trs)
    
    is_source = [0]*len(names)
    is_inter = [0]*len(names)
    is_ecoli = [0]*len(names)
    
    for i in range(len(names)):
        if nodes[i] == source:
            is_source[i] = 1
        elif list(G.predecessors(nodes[i])) == []:
            is_inter[i] = 1
    
    path_formated = pd.DataFrame()
    path_formated['Path ID'] = path_id
    path_formated['Unique ID'] = unique_id
    path_formated['Rule ID'] = rule_id
    path_formated['Left'] = left
    path_formated['Right'] = right
    
    cmpd_formated = pd.DataFrame()
    cmpd_formated['Compound ID']= names
    cmpd_formated['Structure'] = nodes
    cmpd_formated['Source'] = is_source
    cmpd_formated['Intermediate'] = is_inter
    
    return path_formated, cmpd_formated

def clean_path(path_mc):
    n_paths_mc = len(set(path_mc['Path ID']))
    
    p = list(path_mc.groupby(path_mc['Path ID']))
    for i in range(len(p)):
        p[i][1].reset_index(inplace=True)
        p[i][1].drop(labels=['Path ID', 'index'],inplace=True, axis=1)
        p[i] = p[i][1]
    
    to_check = list(range(len(p)))
    repeated = []
    for i in range(len(p)):
        if i not in repeated:
            for j in to_check:
                if i != j and len(p[i]) == len(p[j]):
                    a = p[i].compare(p[j])
                    if len(a) == 0:
                        repeated.append(j)
            for r in repeated:
                try:
                    to_check.remove(r)
                except:
                    continue
    repeated.sort()
    
    new_path_mc = path_mc[path_mc['Path ID'] == 0]
    for i in range(1, n_paths_mc):
        if i not in repeated:
            new_path_mc = pd.concat([new_path_mc, path_mc[path_mc['Path ID'] == i]])
    new_path_mc.reset_index(inplace=True)
    new_path_mc.drop('index', inplace=True, axis=1)
    old_id = new_path_mc['Path ID'].tolist()
    old_id = list(dict.fromkeys(old_id))
    new_id = []
    for i in range(len(old_id)):
        for j in range(len(new_path_mc)):
            if new_path_mc['Path ID'][j] == old_id[i]:
                new_id.append(i)            
    new_path_mc['Path ID'] = new_id
    
    return new_path_mc

pair_name = 'D0P27'
os.makedirs('MC', exist_ok=True)
os.makedirs('MC/'+pair_name, exist_ok=True)

res_formated = scope_format(pair_name)
res_formated.to_csv('MC/'+pair_name+'/MC_Scope.csv', index=False, quoting=csv.QUOTE_ALL, lineterminator='\n')
path_formated, cmpd_formated = path_cmpd_format(pair_name)
path_formated = clean_path(path_formated)
cmpd_formated.to_csv('MC/'+pair_name+'/MC_Compounds.tsv', index=False, sep='\t', lineterminator='\n')
path_formated.to_csv('MC/'+pair_name+'/MC_Pathways.csv', index=False, quoting=csv.QUOTE_ALL, lineterminator='\n')