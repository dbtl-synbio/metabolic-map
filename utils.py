import pandas as pd

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors

def mol_normalizer(smile, mode='inchi'):
    RDLogger.DisableLog('rdApp.*')
    m = Chem.MolFromSmiles(smile)
    Chem.SanitizeMol(m)
    inchi = Chem.MolToInchi(m)
    if mode == 'smile':
        m = Chem.MolFromInchi(inchi)
        smile = Chem.CanonSmiles(Chem.MolToSmiles(m))
        return smile
    return inchi

def mol_similarity(m1, m2, mode='inchi'):
    RDLogger.DisableLog('rdApp.*')
    fpgen = AllChem.GetRDKitFPGenerator()
    if mode == 'smile':
        m1 = Chem.MolFromSmiles(m1)
        m2 = Chem.MolFromSmiles(m2)
    elif mode == 'inchi':
        m1 = Chem.MolFromInchi(m1)
        m2 = Chem.MolFromInchi(m2)
    inchi1 = Chem.MolToInchi(m1).split('/')
    inchi2 = Chem.MolToInchi(m2).split('/')
    c_in1 = ''
    c_in2 = ''
    length = 4
    if len(inchi1) < length or len(inchi2) < length:
        length = min(len(inchi1), len(inchi2))
    for i in range(length):
        c_in1 += inchi1[i]
        c_in2 += inchi2[i]
    if c_in1 == c_in2 and c_in1 != '':
        return 1
    m1_wt = Descriptors.MolWt(Chem.RemoveHs(m1))
    m2_wt = Descriptors.MolWt(Chem.RemoveHs(m2))
    score = DataStructs.cDataStructs.TanimotoSimilarity(fpgen.GetFingerprint(m1), fpgen.GetFingerprint(m2))
    if score == 1:
        score = score - abs(m1_wt-m2_wt)/(m1_wt+m2_wt)
    return score

def smile_in_sink(smile, sink):
    if smile in set(sink):
        return True
    else:
        for s in sink:
            if mol_similarity(smile, s, mode='smile') == 1:
                return True
    return False

def list_smile_in_sink(smiles, sink):
    RDLogger.DisableLog('rdApp.*')
    fpgen = AllChem.GetRDKitFPGenerator()
    f_sink = []
    for i in sink:
        try:
            m = Chem.MolFromSmiles(i)
            f_sink.append(fpgen.GetFingerprint(m))
        except:
            continue
    f_smiles = []
    for i in smiles:
        try:
            m = Chem.MolFromSmiles(i)
            f_smiles.append(fpgen.GetFingerprint(m))
        except:
            continue
    in_sink = [0]*len(smiles)
    for i in range(len(f_smiles)):
        for f in f_sink:
            if DataStructs.cDataStructs.TanimotoSimilarity(f_smiles[i], f) == 1:
                in_sink[i] = 1
                break
    return in_sink

def rule_id_converter(rule, converter, diameter=16):
    r = pd.read_csv('/home/hector/Descargas/retrorules_rr02_rp2_hs/retrorules_rr02_rp2_flat_all.csv')
    rr_id = []
    mnx_id = []
    for i in range(len(r)):
        if r['Diameter'][i] == diameter:
            rr_id.append(r['Rule ID'][i])
            mnx_id.append(r['Legacy ID'][i])
    if type(rule) is list:
        converted_rule = []
        if converter == 'rr':
            for r in rule:
                c = []
                idx = [i for i, x in enumerate(mnx_id) if x == r]
                for i in idx:
                    c.append(rr_id[i])
                converted_rule.append(c)
        elif converter == 'mnx':
            for r in rule:
                c = []
                idx = [i for i, x in enumerate(rr_id) if x == r]
                for i in idx:
                    c.append(mnx_id[i])
                converted_rule.append(c)
        return converted_rule
    else:
        converted_rule = []
        if converter == 'rr':
            idx = [i for i, x in enumerate(mnx_id) if x == rule]
            for i in idx:
                converted_rule.append(rr_id[i])
        elif converter == 'mnx':
            idx = [i for i, x in enumerate(rr_id) if x == rule]
            for i in idx:
                converted_rule.append(mnx_id[i])
        return converted_rule