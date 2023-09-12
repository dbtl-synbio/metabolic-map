import pandas as pd
import networkx as nx

def get_bif(f_ruta, f_scope):
    bif_path = []
    df_r = pd.read_csv(f_ruta)
    df = pd.read_csv(f_scope, usecols =['Substrate','Product','Usage'])
    ini = df_r['A'][0]
    fin = df_r[df_r['Pathway'][0][-1]][0]
    mid =  list(get_scope(len(df_r.iloc[-1]), f_ruta))
    mid.remove(ini)
    prod = df['Product'].tolist()
    padre = df['Substrate'].tolist()
    usage = df['Usage'].tolist()
    if 'retro' in usage:
        for j in range(len(prod)):
            if usage[j] != 'forward' and padre[j] in mid and '.' in prod[j]:
                p_temp = prod[j].split('.')
                if padre[j] not in p_temp and fin not in p_temp:
                    bif_path.append(f_ruta+'$'+padre[j]+'$'+prod[j])
    bif_path = list(dict.fromkeys(bif_path))
    for i in range(len(bif_path)):
        bif_path[i] = bif_path[i].split('$')
    return bif_path
  
def expand_path(G, b_path, ini, filepath):
    alphabet = list(map(chr, range(65,91)))
    df = pd.read_csv(filepath+b_path[0])
    keys = list(df.keys()[1:])
    df = df[keys]
    ruta = []
    for i in keys:
        ruta.append(df[i].tolist())
    b = b_path[2].split('.')
    idx_path = []
    for j in b:
        for i in range(len(df)):
            for k in range(len(ruta)-1):
                if ruta[k][i] == j and ruta[k+1][i] == b_path[1] and '(' not in ruta[0][i]:
                    idx_path.append([i,keys[k],keys[k+1]])
    fin = []
    met = []
    for i in idx_path:
        for j in range(len(ruta)):
            met.append(ruta[j][i[0]])
    met = list(dict.fromkeys(met))
    for i in b:
        if i not in met:
            fin.append(i)
            op = 0
    if len(fin) == 0:
        fin = b
        op = 1
    paths = []
    for i in range(len(ini)):
        for j in range(len(fin)):
            if nx.is_simple_path(G, [ini[i],fin[j]]):
                a = nx.all_simple_paths(G, ini[i], fin[j], cutoff=1)
                for k in a:
                    paths.append(k)
    paths.sort(key = len)
    for i in paths:
        for j in idx_path:
            temp_path = list(df.iloc[j[0]])
            if i[1] in temp_path and op == 1:
                continue
            for t in range(len(temp_path)):
                temp_path[t] = str(temp_path[t])
            t = 0
            if temp_path.count('nan') != 0:
                t = temp_path.index('nan')
            if t != 0:
                p_letra = alphabet[t-1]
                h_letra = alphabet[t]
            else:
                p_letra = alphabet[len(temp_path)-1]
                h_letra = alphabet[len(temp_path)]
            if p_letra not in keys:
                ruta.append([float('nan')]*len(ruta[0]))
                ruta.append([float('nan')]*len(ruta[0]))
                keys.append(p_letra)
                keys.append(h_letra)
            elif h_letra not in keys:
                ruta.append([float('nan')]*len(ruta[0]))
                keys.append(h_letra)
            idx_letra = ruta[0][j[0]].index(j[1])
            str_ruta = ruta[0][j[0]][0:idx_letra+1] +'+('+p_letra+'>'+h_letra+')'+ruta[0][j[0]][idx_letra+1:]
            ruta[0].append(str_ruta)
            for k in range(len(ruta)-1):
                if keys[k+1] == p_letra:
                    ruta[k+1].append(i[0])
                elif keys[k+1] == h_letra:
                    ruta[k+1].append(i[1])
                else:
                    ruta[k+1].append(ruta[k+1][j[0]])
    df = pd.DataFrame(ruta[0],columns=keys[0:1])
    for k in range(len(keys)-1):
        df[keys[k+1]] = ruta[k+1]
    return df

def combine_paths(f):
    df = pd.read_csv(f)
    keys = list(df.keys()[1:])
    df = df[keys]
    alphabet = list(map(chr, range(65,91)))
    ruta = []
    for i in keys:
        ruta.append(df[i].tolist())
    all_count = []
    for i in range(len(df)):
        count = [i]
        if '+' not in ruta[0][i]:
            continue
        for j in range(len(df)):
            if ruta[0][j].count('>') == ruta[0][i].count('>') and ruta[0][j] != ruta[0][i] and list(df.iloc[j][1:-2]) == list(df.iloc[i][1:-2]):
                count.append(j)
        r_comb = []
        for j in count:
            r_comb.append(ruta[0][j])
        r_type = Counter(r_comb)
        dif = []
        for k in list(r_type.keys()):
            temp = []
            for j in count:
                if ruta[0][j] == k:
                    temp.append(j)
            dif.append(temp)
        if len(list(product(*dif))) != 1:
            all_count.append(list(product(*dif)))
    for i in all_count:
        for a in i:
            new_name, p = path_name_maker(ruta, a)
            ruta_test = [new_name]
            n_reac = ruta[0][a[0]].count('>')
            for j in range(n_reac):
                ruta_test.append(ruta[j+1][a[0]])
            for j in a:
                ruta_test.append(ruta[n_reac+1][j])
                ruta_test.append(ruta[n_reac+2][j])
            while len(ruta_test) < len(ruta):
                ruta_test.append(float('nan'))
            for j in range(p-len(ruta)):
                ruta.append([float('nan')]*len(ruta[0]))
            for j in range(len(ruta_test)):
                ruta[j].append(ruta_test[j])
    for i in range(len(ruta)-1):
        if alphabet[i] not in keys:
            keys.append(alphabet[i])
    df = pd.DataFrame(ruta[0],columns=keys[0:1])
    for k in range(len(keys)-1):
        df[keys[k+1]] = ruta[k+1]
    return df

def path_name_maker(ruta, idx):
    n_reac = ruta[0][idx[0]].count('>')
    alphabet = list(map(chr, range(65,91)))
    new_name = ''
    for i in range(n_reac-1):
        new_name += alphabet[i] + '>'
    new_name += alphabet[n_reac-1]
    p = n_reac
    for i in idx:
        r = ''
        for j in ruta[0][i]:
            if j == '+' or j == '>':
                r += j
        place = r.index('+')
        if place == 0:
            new_name = new_name[0:1] + '+('+alphabet[p] +'>'+alphabet[p+1]+')'+new_name[1:]
        elif place == 1:
            idx_b = new_name.index('B') + 1
            new_name = new_name[0:idx_b] + '+('+alphabet[p] +'>'+alphabet[p+1]+')'+new_name[idx_b:]
        elif place == 2:
            idx_c = new_name.index('C') + 1
            new_name = new_name[0:idx_c] + '+('+alphabet[p] +'>'+alphabet[p+1]+')'+new_name[idx_c:]
        elif place == 3:
            idx_d = new_name.index('D') + 1
            new_name = new_name[0:idx_d] + '+('+alphabet[p] +'>'+alphabet[p+1]+')'+new_name[idx_d:]
        elif place == 4:
            idx_e = new_name.index('E') + 1
            new_name = new_name[0:idx_e] + '+('+alphabet[p] +'>'+alphabet[p+1]+')'+new_name[idx_d:]
        p+= 2
    return new_name, p+1
