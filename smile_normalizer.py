from rdkit import Chem
from Filtro.standardizer import Standardizer

"This function normalizes all the SMILE in the same Canon SMILE format."
def smile_normalizer(smile):
    a = Standardizer()
    mol = list()
    for i in range(len(smile)):
        m = Chem.MolFromSmiles(smile[i])
        Chem.SanitizeMol(m)
        m = a.sequence_rr_legacy(m)
        mol.append(m)
        smile[i] = Chem.CanonSmiles(Chem.MolToSmiles(m))
    return smile, mol
