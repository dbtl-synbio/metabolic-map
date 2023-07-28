from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
from Filtro.standardizer import Standardizer

"This function generates the new list of metabolites produced by the initial data."
"You need to provide the initial Smiles and their steps"
def generator(r,smile, step, leg_id, ec_num):
    RDLogger.DisableLog('rdApp.*')
    productos = list()
    resultado = list()
    smile, mol = smile_normalizer(smile)
    for mm in range(len(mol)):
        for i in range(len(r)):
            try:
                pr = r[i].RunReactants( (mol[mm] ,) )
                if len(pr) > 0:
                    m= pr[0][0]
                    Chem.SanitizeMol(m)
                    if Descriptors.MolWt(m) < 1000:
                        b = Chem.CanonSmiles(Chem.MolToSmiles(m))
                        padre = Chem.CanonSmiles(Chem.MolToSmiles(mol[mm]))
                        c = b.split('.')
                        for k in c:
                            productos.append([k,step[mm]])
                        resultado.append(b+"$"+padre+"$"+leg_id[i]+"$"+ec_num[i])
            except:
                continue
    new_smile = list()
    new_step = list()
    for k in productos:
        if k[0] not in smile:
            a, b = smile_normalizer([k[0]])
            new_smile.append(a)
            new_step.append(k[1])
    return new_smile, new_step
