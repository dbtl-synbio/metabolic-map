import pandas as pd
import networkx as nx
from os import listdir

def create_graph(G, all_prod, file_path):
  a = listdir(file_path)
  for f in a:
    df = pd.read_csv(filepath+f, usecols = ['Producto','Padre','Usage'])
    prod = df['Producto'].tolist()
    padre = df['Padre'].tolist()
    usage = df['Usage'].tolist()
    del df
    for i in range(len(prod)):
      b = str(prod[i]).split('.')
      for j in b:
        if j in all_prod and padre[i] in all_prod and usage[i] == 'forward':
          G.add_edge(padre[i], j)
        elif j in all_prod and padre[i] in all_prod and usage[i] == 'both':
          G.add_edge(padre[i],j)
          G.add_edge(j, padre[i])
        elif j in all_prod and padre[i] in all_prod and usage[i] == 'retro':
          G.add_edge(j, padre[i])
  return G
