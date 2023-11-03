from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

"This function normalizes all the SMILE in the same format. By default the format will be InChI, but it can be changed to SMILE"

def mol_normalizer(smile, mode='inchi'):
    m = Chem.MolFromSmiles(smile)
    Chem.SanitizeMol(m)
    inchi = Chem.MolToInchi(m)
    if mode == 'smile':
        m = Chem.MolFromInchi(inchi)
        smile = Chem.CanonSmiles(Chem.MolToSmiles(m))
        return smile
    return inchi

def mol_similarity(m1, m2, mode='inchi'):
    fpgen = AllChem.GetRDKitFPGenerator()
    if mode == 'smile':
        m1 = Chem.MolFromSmiles(m1)
        m2 = Chem.MolFromSmiles(m2)
    elif mode == 'inchi':
        m1 = Chem.MolFromInchi(m1)
        m2 = Chem.MolFromSInchi(m2)
    score = DataStructs.cDataStructs.TanimotoSimilarity(fpgen.GetFingerprint(m1), fpgen.GetFingerprint(m2))
    return score
