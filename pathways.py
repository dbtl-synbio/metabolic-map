from os import listdir, makedirs
import networkx as nx
import pandas as pd
import json

inicial = pd.read_csv('data/Detectable.csv',usecols=['SMILES'])
inicial = inicial['SMILES'].tolist()

producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
producible = list(dict.fromkeys(producible['SMILES']))

folder = listdir('generator_results')
makedirs('generator_results/Updated_graph', exist_ok=True)
folder.remove('Updated_graph')

paths_done = listdir('generator_results/Updated_graph')
paths_done = set(paths_done)
pathways = []
failed_paths = 0
ongoing_paths = 0
error_paths = 0

for i in folder:
    a = listdir('generator_results/'+i)
    inicial_idx = int(i.replace('D',''))
    for j in a:
        if '.json' in j:
            producible_idx = j.replace('.json','')
            producible_idx = producible_idx[producible_idx.index('P')+1:]
            try:
                with open('generator_results/'+i+'/'+j) as f:
                    data = json.load(f)
                G = nx.node_link_graph(data)
                if nx.get_node_attributes(G, 'score')[inicial[inicial_idx]] == 1:
                    pathways.append(j)
                elif nx.get_node_attributes(G, 'score')[inicial[inicial_idx]] == 0:
                    failed_paths += 1
                else:
                    ongoing_paths += 1
            except:
                error_paths += 1
                continue

print('Pairs done:', len(pathways))
print('Pairs with no solution:', failed_paths)
print('Pairs ongoing:', ongoing_paths)
print('Errors:', error_paths)
df = pd.DataFrame(pathways, columns=['Name'])
df.to_csv('data/Pathways.csv', mode='a', index=False)