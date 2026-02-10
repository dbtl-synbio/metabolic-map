import pandas as pd
import networkx as nx
import json
from os import listdir
from utils import mol_similarity

def graph_updater(res_filename, graph_filename, save_filename, target):
    df = pd.read_csv(res_filename)
    resultado=df['Resultado'].to_list()
    resultado=list(dict.fromkeys(resultado))
    try:
        resultado.remove('Resultado')
    except:
        a = 1
    padre = []
    reac = []
    prod = []
    ec_num = []
    for i in resultado:
        a,b,c,d= i.split('$')
        padre.append(b)
        reac.append(c)
        prod.append(a)
        ec_num.append(d)
    del(a,b,c,d,df,i,resultado)
    
    data = json.load(open(graph_filename))
    G = nx.node_link_graph(data)
    a = nx.get_node_attributes(G, 'score')
    for i in a:
        if a[i] == 1:
            if mol_similarity(target, i, mode='smile') == 1:
                producible = i
        if a[i] == 1 and list(G.predecessors(i)) == []:
            inicial = i
            
    for i in nx.shortest_simple_paths(G, inicial, producible):
        path = i
        
    rules_id = []
    for i in range(len(path)-1):
        for j in range(len(prod)):
            if padre[j] == path[i] and prod[j] == path[i+1]:
                rules_id.append([padre[j], prod[j], reac[j]])

    other = []
    for i in range(len(rules_id)):
        for j in range(len(prod)):
            if rules_id[i][0] == padre[j] and rules_id[i][2] == reac[j] and prod[j] not in path:
                other.append([padre[j], prod[j], reac[j]])
                
    H = nx.DiGraph()
    graph_info = rules_id + other
    for i in graph_info:
        H.add_node(i[0], type='chemical')
        H.add_node(i[1], type='chemical')
        H.add_node(i[2], type='reaction')
        H.add_edge(i[2], i[0])
        H.add_edge(i[1], i[2])
    data = nx.node_link_data(H)
    with open(save_filename, 'w') as f:
        f.write(json.dumps(data))

list_graph = pd.read_csv('data/Pathways.csv')
list_graph = list_graph[list_graph['Name'] != 'Name']
list_graph.reset_index(inplace=True, drop=True)
list_graph.to_csv('data/Pathways.csv')
list_graph = list_graph['Name'].tolist()
done_graph = set(listdir('generator_results/Updated_graph'))

producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
producible = producible['SMILES'].tolist()

todo = [x for x in list_graph if x not in done_graph]
for i in todo:
    name = i.split('_')[1]
    name = name.replace('.json','')
    res_filename = 'generator_results/'+name.split('P')[0]+'/Results_'+name+'.csv'
    graph_filename = 'generator_results/'+name.split('P')[0]+'/Graph_'+name+'.json'
    save_filename = 'generator_results/Updated_graph/Graph_'+name+'.json'
    target = producible[int(name.split('P')[1])]
    graph_updater(res_filename, graph_filename, save_filename, target)
