import pandas as pd
import numpy as np
import os
from rdkit import Chem, RDLogger
import csv

from utils import list_smile_in_sink, mol_similarity

def merger(pair_name):
    # Leemos el archivo MC
    path_mc = pd.read_csv('MC/'+pair_name+'/MC_Pathways.csv')
    cmpd_mc = pd.read_csv('MC/'+pair_name+'/MC_Compounds.tsv',sep='\t')
    scope_mc = pd.read_csv('MC/'+pair_name+'/MC_Scope.csv')
    
    right_mc = []
    for i in range(len(path_mc)):
        idx = path_mc['Right'][i].index('.')
        right_mc.append(path_mc['Right'][i][idx+1:])
    right_mc = list(dict.fromkeys(right_mc))
    # Encontramos los intermedios
    inter = pd.read_csv('data/Intermediate.csv')
    inter= inter['Smile'].tolist()
    idx = []
    for i in range(len(cmpd_mc)):
        if cmpd_mc['Structure'][i] in inter and cmpd_mc['Compound ID'][i] not in right_mc:
            idx.append(inter.index(cmpd_mc['Structure'][i]))
    # Preparamos los archivos ME
    scope_me = pd.concat([scope_mc])
    path_me = pd.concat([path_mc])
    cmpd_me = pd.concat([cmpd_mc])
    
    inter_targets = []
    
    for k in idx:
        # Leemos los archivos intermedios
        try:
            path_in = pd.read_csv('Intermedia_path/source_I'+str(k)+'/outpath/out_paths.csv')
            cmpd_in = pd.read_csv('Intermedia_path/source_I'+str(k)+'/outpath/compounds.txt',sep='\t')
            cmpd_in.drop_duplicates(ignore_index=True, inplace=True)
            scope_in = pd.read_csv('Intermedia_path/source_I'+str(k)+'/Inter_'+str(k)+'_scope.csv')
        except:
            continue
        
        # Corregimos las SMILES de Compounds_IN
        for i in range(len(cmpd_in)):
            cmpd_in['Structure'][i] = Chem.CanonSmiles(cmpd_in['Structure'][i])
    
        # Obtenemos la equivalencia de Compounds ID
        c_id_in = cmpd_in['Compound ID'].tolist()
        c_id_mc = cmpd_mc['Compound ID'].tolist()
        s_id_in = cmpd_in['Structure'].tolist()
        s_id_mc = cmpd_mc['Structure'].tolist()
    
        c_equivalence = []
        cmpd_idx = len(s_id_mc)+1
        for i in range(len(s_id_in)):
            try:
                idx_eq = s_id_mc.index(s_id_in[i])
                c_equivalence.append(c_id_mc[idx_eq])
            except:
                idx_eq = -1
                for j in range(len(s_id_mc)):
                    if mol_similarity(s_id_mc[j], s_id_in[i], mode='smile') == 1:
                        idx_eq = j
                if idx_eq == -1:
                    if 'CMPD' in c_id_in[i]:
                        c_equivalence.append('CMPD_'+str(cmpd_idx).zfill(10))
                        cmpd_idx += 1
                    else:
                        c_equivalence.append(c_id_in[i])
                else:
                    c_equivalence.append(c_id_mc[idx_eq])
                    
        # Corregimos los Pathways_IN para evitar que el TARGET sea parte del pathway
        checked_path_in = pd.DataFrame({'Path ID':[],'Unique ID':[], 'Rule ID':[], 'Left':[], 'Right':[]})
    
        for i in range(len(path_in)):
            left = path_in['Left'][i].split(':')
            l = []
            for j in range(len(left)):
                idx = c_id_in.index(left[j][2:])
                l.append(c_equivalence[idx])
            idx = c_id_in.index(path_in['Right'][i][2:])
            right = path_in['Right'][i][:2] + c_equivalence[idx]
            new_left = ''
            for j in range(len(l)):
                new_left += left[j][:2] + l[j] +':'
            if new_left[-1] == ':':
                new_left = new_left[:-1]
            if 'TARGET_0000000001' in l:
                continue
            new_row = pd.DataFrame({'Path ID':path_in['Path ID'][i],'Unique ID':path_in['Unique ID'][i].replace('TRS_0','TRS_'+str(k)),'Rule ID':path_in['Rule ID'][i], 'Left':new_left, 'Right':right}, index=[0])
            checked_path_in = pd.concat([checked_path_in, new_row], ignore_index=True)
            
        if len(checked_path_in) == 0:
            continue
        # Corregimos el Scope_IN y lo añadimos al ME
        s = []
        for i in range(len(scope_in)):
            s.append(scope_in['Transformation ID'][i].replace('TRS_0','TRS_'+str(k)))
        scope_in['Transformation ID'] = s
        scope_me = pd.concat([scope_me, scope_in], ignore_index=True)
        
             
        # Calculamos las combinaciones de los pathways
        n_paths_in = len(set(checked_path_in['Path ID']))
        n_paths_mc = len(set(path_me['Path ID']))
        target = c_equivalence[c_id_in.index('TARGET_0000000001')]
        inter_targets.append(target)
        
        duplicate_paths = []
        for i in range(n_paths_mc):
            for j in range(len(path_me)):
                if path_me['Path ID'][j] == i and target in path_me['Left'][j]:
                    duplicate_paths.append(i)
                    break
        unique_paths = []
        duplicate_paths_in = checked_path_in['Path ID'].tolist()
        duplicate_paths_in = list(dict.fromkeys(duplicate_paths_in))
        for i in duplicate_paths_in:
            for j in duplicate_paths:
                unique_paths.append([i,j])
        
        # Obtenemos el nuevo Pathways ME
        new_path_me = pd.DataFrame()
        for i in range(len(unique_paths)):
            sub_path_me = path_me[path_me['Path ID'] == unique_paths[i][1]]
            sub_path_me['Path ID'] = [i] * len(sub_path_me)
            sub_checked_path_in =  checked_path_in[checked_path_in['Path ID'] == unique_paths[i][0]]
            sub_checked_path_in['Path ID'] = [i] * len(sub_checked_path_in)
            new_path_me = pd.concat([new_path_me, sub_path_me, sub_checked_path_in], ignore_index=True)
        path_me = pd.concat([path_me, new_path_me], ignore_index=True)
        path_me.drop_duplicates(inplace=True, ignore_index=True)
                    
        # Obtenemso el nuevo Compound ME
        compound_id = []
        structure = []
        for i in range(len(c_equivalence)):
            for j in range(len(path_me)):
                if c_equivalence[i] in path_me['Left'][j] or c_equivalence[i] in path_me['Right'][j]:
                    compound_id.append(c_equivalence[i])
                    structure.append(s_id_in[i])
                    break
        checked_cmpd_in = pd.DataFrame({'Compound ID': compound_id,'Structure': structure})
        cmpd_me = pd.concat([cmpd_me, checked_cmpd_in], ignore_index=True)
        cmpd_me.drop_duplicates(ignore_index=True, inplace=True)
    # Añadimos la columna Intermediate a Compound ME
    intermediate = [0] * len(cmpd_me)
    for i in range(len(cmpd_me)):
        if cmpd_me['Compound ID'][i] in inter_targets:
            intermediate[i] = 1
    cmpd_me['Intermediate'] = intermediate
    
    return scope_me, path_me, cmpd_me

def add_chassis(chassis_folder, cmpd_me):
    in_sink = []*len(cmpd_me)
    chasis_names = os.listdir(chassis_folder)
    in_chasis = []
    for i in chasis_names:
        chasis = pd.read_csv(chassis_folder+i, usecols=['SMILES'])
        chasis = chasis['SMILES'].tolist()
        chasis = list(dict.fromkeys(chasis))
        good_chasis = []
        for j in chasis:
            if not pd.isna(j):
                good_chasis.append(j)
        in_chasis.append(list_smile_in_sink(cmpd_me['Structure'], good_chasis))
    for i in in_chasis[0]:
        if i == 1:
            in_sink.append([chasis_names[0].replace('.csv','')])
        else:
            in_sink.append([])
    for c in range(1, len(chasis_names)):
        for i in range(len(in_chasis[c])):
            if in_chasis[c][i] == 1:
                in_sink[i].append(chasis_names[c].replace('.csv',''))
                
    cmpd_me['Chassis'] = in_sink
    cmpd_me.drop_duplicates(subset=cmpd_me.columns.difference(['Chassis']), ignore_index=True, inplace=True)
    
    return cmpd_me

RDLogger.DisableLog('rdApp.*')
pair_name = 'D0P27'

scope_me, path_me, cmpd_me = merger(pair_name)

producible = pd.read_csv('data/Producible.csv',usecols=['SMILES'])
producible = list(dict.fromkeys(producible['SMILES']))

cmpd_me = add_chassis('data/chassis/', cmpd_me)

max = 0
for i in range(len(cmpd_me)):
    idx = int(cmpd_me['Compound ID'][i].split('_')[1])
    if 'CMPD' in cmpd_me['Compound ID'][i] and idx > max:
        max = idx
for i in range(len(cmpd_me)):
    if 'EColi' in cmpd_me['Compound ID'][i]:
        old_name = cmpd_me['Compound ID'][i]
        cmpd_me['Compound ID'][i] = 'CMPD_' + str(max+1).zfill(10)
        for j in range(len(path_me)):
            left = path_me['Left'][j].split(':')
            for k in range(len(left)):
                idx = left[k].index('.')
                n_left = left[k][:idx+1]
                if left[k][idx+1:] == old_name:
                    left[k] = n_left + 'CMPD_' + str(max+1).zfill(10)
            path_me['Left'][j] = ''
            for k in left:
                path_me['Left'][j] += k +':'
            path_me['Left'][j] = path_me['Left'][j][:-1]
        max +=1
cmpd = []
for i in range(len(path_me)):
    idx = path_me['Right'][i].index('.') 
    cmpd.append(path_me['Right'][i][idx+1:])
    left = path_me['Left'][i].split(':')
    for k in range(len(left)):
        idx = left[k].index('.')
        left[k] = left[k][idx+1:]
    cmpd += left
    
cmpd = list(dict.fromkeys(cmpd))

remove =[]
for i in range(len(cmpd_me)):
    if cmpd_me['Compound ID'][i] not in cmpd:
        remove.append(i)
cmpd_me.drop(remove, inplace=True)
path_me.sort_values(by=['Path ID','Unique ID'], ignore_index=True, inplace=True)
os.makedirs('ME',exist_ok=True)
os.makedirs('ME/'+pair_name, exist_ok=True)
scope_me.to_csv('ME/'+pair_name+'/ME_Scope.csv', index=False, quoting=csv.QUOTE_ALL, lineterminator='\n')
path_me.to_csv('ME/'+pair_name+'/ME_Pathways.csv', index=False, quoting=csv.QUOTE_ALL, lineterminator='\n')
cmpd_me.to_csv('ME/'+pair_name+'/ME_Compounds.tsv', index=False, sep='\t', lineterminator='\n')