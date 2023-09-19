import pandas as pd

def path_cleaner(filename):
    df = pd.read_csv(filename)
    keys = list(df.keys()[1:])
    df = df[keys]
    alphabet = list(map(chr, range(65,91)))
    ruta = []
    for i in range(len(df)):
        aa = list(df.iloc[i][1:])
        for j in range(len(aa)):
            aa[j] = str(aa[j])
        aa.sort()
        ruta.append(aa)
    duplicate = []
    for i in range(len(ruta)):
        for j in range(len(ruta)):
            r1 = ''
            r2 = ''
            for r in df['Pathway'][i]:
                if r not in alphabet:
                    r1 += r
            for r in df['Pathway'][j]:
                if r not in alphabet:
                    r2 += r
            if ruta[i] == ruta[j] and i != j and r1 == r2:
                duplicate.append([i, j])
    for d in duplicate:
        if [d[1], d[0]] in duplicate:
            duplicate.remove([d[1], d[0]])
    remove = []
    for d in duplicate:
        remove.append(d[1])
    remove = list(dict.fromkeys(remove))
    remove.sort()
    for r in remove:
        df = df.drop(r, axis=0)
    df.to_csv(filename)
