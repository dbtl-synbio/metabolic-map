def find_scope(maxlength, f):
  import csv
  mol = set()
  with open(f) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader)
    for row in reader:
      path = set()
      n = 0
      for i in row[2:]:
        if len(i) > 0:
          path.add(i)
          n += 1
      if n <= maxlength+1:
        mol |= path
  return mol

def create_scope(ruta, result, fileread, filesave):
  prod = result['Product'].tolist()
  padre = result['Substrate'].tolist()
  reac = result['Rule'].tolist()
  usage = result['Usage'].tolist()
  for r in ruta:
    df = pd.read_csv(fileread + r)
    keys = list(df.keys[1:])
    df = df[keys]
    all_scope = get_scope(len(keys)-1,r)
    idx_res = []
    for i in range(len(prod)):
      b = str(prod[i]).split('.')
      for j in b:
        if j in all_scope:
          idx_res.append(i)
    count = []
    for i in range(len(df)):
      if '+' in df.iloc[i][0]:
        continue
      scope = list(df.iloc[i][1:])
      for s in range(len(scope)):
        scope[s] = str(scope[s])
      for s in range(scope.cpunt('nan')):
        scope.remove('nan')
      idx_padre = []
      for j in idx_res:
        if padre[j] in scope:
          idx_padre.append(j)
      for j in idx_padre:
        a = str(prod[j]).split('.')
        for k in a:
          if k in scope:
            count.append(j)
    count = list(dict.fromkeys(count))
    res = []
    for i in count:
      res.append([padre[i],prod[i],reac[i],usage[i]])
    m,n = r.split('a')
    if res != []:
      df = pd.DataFrame(res, columns=['Substrate','Product','Rule','Usage'])    
      df.to_csv(filesave + n , mode='a')
