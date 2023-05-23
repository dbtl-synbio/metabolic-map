import pandas as pd

def get_all_participants(ini_met, f):
  df = pd.read_csv(f, usecols=['Producto','Padre'])
  prod = df['Producto'].tolist()
  padre = df['Padre'].tolist()
  del df
  p = set(ini_met)
  count = []
  for i in range(len(prod)):
    b = str(prod[i]).split('.')
    for j in b:
      if j in p:
        count.append(padre[i])
        break
  count = list(dict.fromkeys(count))
  count2 = []
  for i in count:
    if i not in p:
      count2.append(i)
  ini_met += count2
  while (len(count2) != 0):
    p = set(ini_met)
    count = []
    for i in in range(len(prod)):
      b = str(prod[i]).split('.')
      for j in b:
        if j in p:
          count.append(padre[i])
          break
    count = list(dict.fromkeys(count))
    count2 = []
    for i in count:
      if i not in p:
        count2.append(i)
    ini_met += count2
  return ini_met
