import pandas as pd
import networkx as nx

def find_path(G, ini, fin):
  path = []
  for i in range(len(ini)):
    for j in range(len(fin)):
      if G.has_node(ini[i]) and G.has_node(fin[j]):
        if nx.has_path(G, ini[i], fin[j]):
          path.append([i,j])
  return path

def generate_path(G, ini, fin, limit=100):
  list_paths = []
  a = nx.shortest_simple_paths(G, ini, fin)
  for i in a:
    list_paths.append(i)
    if len(list_paths) > (limit -1):
      break
  df = pd.DataFrame(list_paths)
  df_len = len(df.iloc[-1])
  keys = list(map(chr, range(65,91)))
  df.columns = keys[0:df_len]
