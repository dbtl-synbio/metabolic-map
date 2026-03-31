import pandas as pd
from os import listdir
import csv
import os

pair_name = 'D0P27'

path_mc = pd.read_csv('MC/'+pair_name+'/MC_Pathways.csv')
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
new_path_mc.to_csv('/home/hector/Doctorado/MonteCarlo/pair_data/'+pair_name+'/MC_Pathways_mod.csv', index=False, quoting=csv.QUOTE_ALL, lineterminator='\n')