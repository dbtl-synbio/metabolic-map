import os
import sys
import json
import argparse

sys.path.append(os.path.join(os.getenv("HOME"),'Doctorado/scripts/rptools'))
sys.path.append(os.path.join(os.getenv("HOME"),'Doctorado/scripts'))

from rptools.rpviz.utils import annotate_cofactors, annotate_chemical_svg, get_autonomous_html, parse_all_pathways
from rptools.rpviz.Viewer import Viewer

def get_one_pathways(path, scope):
    path = path.reset_index(drop=True)
    pathway = {'Path_id':'Path_'+str(path['Path ID'][0]), 'node_ids':[], 'edge_ids':[],'nb_steps':10,'scores':{'global_score':0}}
    
    # Get all node ids
    node_ids = []
    for i in range(len(path)):
        node_ids.append(path['Rule ID'][i])
        node_ids += path['Left'][i].split(':')
        node_ids.append(path['Right'][i])
    for i in range(len(node_ids)):
        try:
            idx = node_ids[i].index('.')
            node_ids[i] = node_ids[i][idx+1:]
        except:
            continue
    node_ids = list(dict.fromkeys(node_ids))
    
    # Get all edges ids
    edge_ids = []
    for i in range(len(path)):
        left = path['Left'][i].split(':')
        for j in left:
            idx = j.index('.')
            edge_ids.append(j[idx+1:]+'_'+path['Rule ID'][i])
        idx = path['Right'][i].index('.')
        edge_ids.append(path['Rule ID'][i]+'_'+path['Right'][i][idx+1:])
    
    max_score = 0
    t_ids = scope['Transformation ID'].tolist()
    for i in path['Unique ID']:
        idx= t_ids.index(i[:-2])
        if scope['Score'][idx] > max_score:
            max_score = scope['Score'][idx]
    pathway['node_ids'] = node_ids
    pathway['edge_ids'] = edge_ids
    pathway['nb_steps'] = len(path)-1
    pathway['scores']['global_score'] = max_score
    return pathway

def get_detectable_info(nodes, cmpd):
    import pandas as pd
    from utils import mol_similarity
    detectable = pd.read_csv('../data/Sensbio_info.csv')
    inicial = pd.read_csv('../data/Detectable.csv',usecols=['Name','SMILES'])
    for i in range(len(cmpd)):
        nodes[i]['data']['species'] = []
        nodes[i]['data']['tf'] = []
        nodes[i]['data']['biblio_ref'] = []
        nodes[i]['data']['data_ref'] = []
        nodes[i]['data']['ncbi'] = []
        nodes[i]['data']['uniprot'] = []
        nodes[i]['data']['aa_seq'] = []
        if 'TARGET' in cmpd['Compound ID'][i]:
            idx = inicial['SMILES'].tolist().index(cmpd['Structure'][i])
            for j in range(len(detectable)):
                if detectable['Molecule'][j] == inicial['Name'][idx]:
                    nodes[i]['data']['species'].append(detectable['Species'][j])
                    nodes[i]['data']['tf'].append(detectable['TF'][j])
                    nodes[i]['data']['biblio_ref'].append(detectable['Bibliographic_ref'][j])
                    nodes[i]['data']['data_ref'].append(detectable['Database_ref'][j])
                    nodes[i]['data']['ncbi'].append(detectable['NCBI_Accession'][j])
                    nodes[i]['data']['uniprot'].append(detectable['UniProt'][j])
                    nodes[i]['data']['aa_seq'].append(detectable['AA_sequence'][j])
    return nodes

def get_producible_info(nodes, cmpd):
    import pandas as pd
    from utils import mol_similarity
    producible = pd.read_csv('../data/Producible.csv')
    for i in range(len(cmpd)):
        nodes[i]['data']['pmid'] = None
        nodes[i]['data']['keggid'] = None
        nodes[i]['data']['mnxid'] = None
        if cmpd['Source'][i] == 1:
            for p in range(len(producible['SMILES'])):
                try:
                    if mol_similarity(cmpd['Structure'][i], producible['SMILES'][p], mode='smile') == 1:
                        nodes[i]['data']['pmid'] = producible['PMID1'][p]
                        nodes[i]['data']['keggid'] = producible['KEGG_ID'][p]
                        nodes[i]['data']['mnxid'] = producible['MNX_ID'][p]
                except:
                    continue
    return nodes

def get_main_compounds(nodes, pathways):
    pathways.drop(labels=['Path ID'], inplace=True, axis=1)
    pathways.drop_duplicates(inplace=True, ignore_index=True)
    cmpd = []
    for i in range(len(pathways)):
        if 'TRS_0' in pathways['Unique ID'][i]:
            idx = pathways['Right'][i].index('.') 
            cmpd.append(pathways['Right'][i][idx+1:])
            left = pathways['Left'][i].split(':')
            for k in range(len(left)):
                idx = left[k].index('.')
                left[k] = left[k][idx+1:]
            cmpd += left
    for i in range(len(nodes)):
        if nodes[i]['data']['id'] in cmpd:
            nodes[i]['data']['main_path'] = True
        else:
            nodes[i]['data']['main_path'] = False
    return nodes
    
def get_nodes(compounds, pathways, scope, prod_path=False):
    from rdkit import Chem, RDLogger
    from utils import rule_id_converter
    RDLogger.DisableLog('rdApp.*')
    nodes = []
    for i in range(len(compounds)):
        nodes.append({'data':{'id':compounds['Compound ID'][i]}})
        nodes[i]['data']['path_ids'] = []
        nodes[i]['data']['inter_ids'] = []
        nodes[i]['data']['prod_path'] = prod_path
        nodes[i]['data']['label'] = compounds['Compound ID'][i]
        nodes[i]['data']['all_labels'] = compounds['Compound ID'][i]
        nodes[i]['data']['type'] = 'chemical'
        nodes[i]['data']['smiles'] = compounds['Structure'][i]
        try:
            m = Chem.MolFromSmiles(compounds['Structure'][i])
            nodes[i]['data']['inchi'] = Chem.MolToInchi(m)
            nodes[i]['data']['inchikey'] = Chem.MolToInchiKey(m)
        except:
            nodes[i]['data']['inchi'] = ''
            nodes[i]['data']['inchikey'] = ''
        if 'TARGET' in compounds['Compound ID'][i]:
            nodes[i]['data']['target_chemical'] = True
            nodes[i]['data']['short_label'] = 'TARGET'
        else:
            nodes[i]['data']['target_chemical'] = False
            if 'CMPD_' in compounds['Compound ID'][i]:
                nodes[i]['data']['short_label'] = 'C' + str(int(compounds['Compound ID'][i].split('_')[1]))
            else:
                nodes[i]['data']['short_label'] = compounds['Compound ID'][i]
        nodes[i]['data']['chassis'] = compounds['Chassis'][i]
        if compounds['Intermediate'][i] == 0:
            nodes[i]['data']['inter_chemical'] = False
        else:
            nodes[i]['data']['inter_chemical'] = True
        if compounds['Source'][i] == 0:
            nodes[i]['data']['source_chemical'] = False
        else:
            nodes[i]['data']['source_chemical'] = True
            nodes[i]['data']['prod_path'] = True
        nodes[i]['data']['cofactor'] = False
        nodes[i]['data']['svg'] = ''
        nodes[i]['data']['rsmiles'] = None
        nodes[i]['data']['rule_ids'] = None
        nodes[i]['data']['rule_score'] = None
        nodes[i]['data']['rxn_template_ids'] = None
        nodes[i]['data']['ec_numbers'] = None
        nodes[i]['data']['xlinks'] = [{'db_name':'','entity_id':compounds['Compound ID'][i],'url':''}]
        nodes[i]['data']['sink_chemical'] = False
    nodes = get_producible_info(nodes, compounds)
    nodes = get_detectable_info(nodes, compounds)
    nodes = get_main_compounds(nodes, pathways)
    for i in range(len(pathways)):
        nodes.append({'data':{'id':pathways['Rule ID'][i]}})
        nodes[-1]['data']['path_ids'] = []
        nodes[-1]['data']['inter_ids'] = []
        nodes[-1]['data']['prod_path'] = prod_path
        nodes[-1]['data']['label'] = pathways['Rule ID'][i]
        nodes[-1]['data']['all_labels'] = pathways['Rule ID'][i]
        nodes[-1]['data']['short_label'] = pathways['Rule ID'][i]
        nodes[-1]['data']['type'] = 'reaction'
        nodes[-1]['data']['inchi'] = None
        nodes[-1]['data']['inchikey'] = None
        nodes[-1]['data']['cofactor'] = None
        for j in range(len(scope)):
            if scope['Transformation ID'][j] == pathways['Unique ID'][i][:-2]:
                nodes[-1]['data']['ec_numbers'] = scope['EC number'][j]
                nodes[-1]['data']['rule_score'] = scope['Score'][j]
                nodes[-1]['data']['rsmiles'] = scope['Reaction SMILES'][j]
        nodes[-1]['data']['rule_ids'] = pathways['Rule ID'][i]
        nodes[-1]['data']['rxn_template_ids'] = pathways['Rule ID'][i]
        nodes[-1]['data']['chassis'] = None
        nodes[-1]['data']['inter_chemical'] = None
        nodes[-1]['data']['source_chemical'] = None
        nodes[-1]['data']['target_chemical'] = None
        nodes[-1]['data']['smiles'] = None
        nodes[-1]['data']['svg'] = None
        nodes[-1]['data']['pmid'] = None
        nodes[-1]['data']['keggid'] = None
        nodes[-1]['data']['mnxid'] = None
           
    return nodes

def get_edges(pathways, prod_path=False):
    edges = []
    for i in range(len(pathways)):
        left = pathways['Left'][i].split(':')
        for j in left:
            idx = j.index('.')
            edges.append({'data':{'id':j[idx+1:]+'_'+pathways['Rule ID'][i]}})
            edges[-1]['data']['path_ids'] = []
            edges[-1]['data']['inter_ids'] = []
            edges[-1]['data']['prod_path'] = prod_path
            edges[-1]['data']['source'] = j[idx+1:]
            edges[-1]['data']['target'] = pathways['Rule ID'][i]
        idx = pathways['Right'][i].index('.')
        edges.append({'data':{'id':pathways['Rule ID'][i] +'_'+pathways['Right'][i][idx+1:]}})
        edges[-1]['data']['path_ids'] = []
        edges[-1]['data']['inter_ids'] = []
        edges[-1]['data']['prod_path'] = prod_path
        edges[-1]['data']['source'] = pathways['Rule ID'][i]
        edges[-1]['data']['target'] = pathways['Right'][i][idx+1:]
        
    return edges

def insert_paths_ids(network, pathways_info):
    for p in pathways_info:
        for i in network['elements']['nodes']:
            if i['data']['id'] in pathways_info[p]['node_ids']:
                i['data']['path_ids'].append(p)
        for i in network['elements']['edges']:
            if i['data']['id'] in pathways_info[p]['edge_ids']:
                i['data']['path_ids'].append(p)
                i['data']['path_ids'] = list(dict.fromkeys(i['data']['path_ids']))
    return network

def get_intermediate_paths(network, pathways):
    for i in range(len(pathways)):
        inter_node = []
        inter_edge = []
        if 'TRS_0' not in pathways['Unique ID'][i]:
            inter_node.append(pathways['Rule ID'][i])
            idx = pathways['Right'][i].index('.')
            inter_node.append(pathways['Right'][i][idx+1:])
            inter_edge.append(pathways['Rule ID'][i] +'_'+pathways['Right'][i][idx+1:])
            left = pathways['Left'][i].split(':')
            for j in left:
                idx = j.index('.')
                inter_node.append(j[idx+1:])
                inter_edge.append(j[idx+1:] +'_'+pathways['Rule ID'][i])
        if inter_node != []:
            inter_id = pathways['Unique ID'][i].split('_')
            for j in network['elements']['nodes']:
                if j['data']['id'] in inter_node:
                    j['data']['inter_ids'].append(inter_id[0]+'_'+inter_id[1])
            for j in network['elements']['edges']:
                if j['data']['id'] in inter_edge:
                    j['data']['inter_ids'].append(inter_id[0]+'_'+inter_id[1])
    return network

def host_pathways(network, compounds):
    import pandas as pd
    from utils import mol_similarity
    
    index_prod = pd.read_csv('../data/Producible.csv', usecols=['SMILES'])
    index_prod = index_prod['SMILES'].tolist()
    for i in range(len(compounds)):
        if compounds['Source'][i] == 1:
            source = compounds['Structure'][i]
            source_name = compounds['Compound ID'][i]
    
    idx_prod = -1
    for i in range(len(index_prod)):
        try:
            if mol_similarity(source, index_prod[i], mode='smile') == 1:
                idx_prod = i
        except:
            continue
    if idx_prod == -1:
        return network
    
    try:
        path = pd.read_csv('../data/Producible_Pathways/P'+str(idx_prod)+'_Pathways.csv')
        cmpd = pd.read_csv('../data/Producible_Pathways/P'+str(idx_prod)+'_Compounds.tsv', sep='\t')
    except:
        return network
    
    cmpd['Source'] = [0]*len(cmpd)
    cmpd['Intermediate'] = [0]*len(cmpd)
    scope = pd.DataFrame()
    scope['Transformation ID'] = path['Unique ID']
    scope['Score'] = [1]*len(scope)
    scope['EC number'] = ['']*len(scope)
    
    nodes = get_nodes(cmpd, path, scope, prod_path = True)
    edges = get_edges(path, prod_path = True)
    
    prod_name = ''
    for i in range(len(cmpd['Structure'])):
        try:
            if mol_similarity(source, cmpd['Structure'][i], mode='smile') == 1:
                prod_name = cmpd['Compound ID'][i]
        except:
            continue
    if prod_name == '':
        return network
    
    nodes_corr = []
    for n in nodes:
        if n['data']['id'] != prod_name:
            nodes_corr.append(n)
    for e in range(len(edges)):
        if edges[e]['data']['source'] == prod_name:
            edges[e]['data']['source'] = source_name
            edges[e]['data']['id'] = source_name +'_'+ edges[e]['data']['target']
        elif edges[e]['data']['target'] == prod_name:
            edges[e]['data']['target'] = source_name
            edges[e]['data']['id'] = edges[e]['data']['source'] +'_'+ source_name
    
    network['elements']['nodes'] += nodes_corr
    network['elements']['edges'] += edges
    
    return network

def get_all_pathways(input_files):
    network = {'elements': {'nodes': [], 'edges': []}}
    pathways_info = {}
    
    import pandas as pd
    for i in input_files:
        if 'Scope.csv' in i:
            scope = pd.read_csv(i)
        elif 'Compounds.tsv' in i:
            compounds = pd.read_csv(i, sep='\t')
        elif 'Pathways.csv' in i:
            pathways = pd.read_csv(i)
        else:
            raise NotImplementedError('Some files on the input folder do not follow the format. Exit')
    
    # Get all pathways
    for i in list(set(pathways['Path ID'])):
        pathways_info['Path_'+str(i)] = get_one_pathways(pathways[pathways['Path ID'] == i], scope)
    # Get nodes
    network['elements']['nodes'] = get_nodes(compounds, pathways, scope)
    # Get edges
    network['elements']['edges'] = get_edges(pathways)
    
    #Add path ids to network
    network = insert_paths_ids(network, pathways_info)
    
    #Add intermediate paths ids to network
    network = get_intermediate_paths(network, pathways)
    network = host_pathways(network, compounds)
    return network, pathways_info

def __run(args):
    input_file = '../ME/' + args
    input_files = [os.path.join(input_file, file) for file in os.listdir(input_file)]
    
    network, pathways_info = get_all_pathways(input_files)
    output_folder = '../json_pair_files/P' + str(args.split('P')[1])
    os.makedirs(output_folder, exist_ok=True)
    
    network_out_file = os.path.join(output_folder, args+'_network.json')
    path_out_file = os.path.join(output_folder, args+'_pathway.json')
    with open(network_out_file,'w') as f:
        f.write(json.dumps(network))
    with open(path_out_file,'w') as f:
        f.write(json.dumps(pathways_info))
        
if __name__ == '__main__':
    args = os.listdir('../ME')
    for i in args:
        det, prod = i.split('P')
        if os.path.isdir('../json_pair_files/P'+prod+'/') == True:
            a = os.listdir('../json_pair_files/P'+prod+'/')
            if det+'P'+prod+'_network.json' in a:
                continue
        #__run(i)
        try:
            __run(i)
            print(i+' completed')
        except:
            print(i +' not working')