import gzip
import json
import typing as tp

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions


def load_data(
        path: str,
) -> tp.Tuple[tp.List, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns the building blocks, the scaffolds, the reactions, and the constant regions.
    """
    scaffolds = pd.read_excel(path, sheet_name='scaffolds')
    reactions = pd.read_excel(path, sheet_name='smarts')
    consts = pd.read_excel(path, sheet_name='const')
    bbs = []
    try:
        while True:
            step = len(bbs) + 1
            bbs += [pd.read_excel(path, sheet_name=f'step{step}').fillna('')]
    except:
        return (bbs, scaffolds, reactions, consts)


def get_smiles(
        scaffold_id: str,
        scaffolds: pd.DataFrame,
) -> str:
    mask = scaffolds['ScaffoldID'] == scaffold_id
    return scaffolds['SMILES'][mask].values[0]


def get_smarts(
        reaction_type: str,
        reactions: pd.DataFrame,
) -> str:
    try:
        mask = reactions['ReactionType'] == reaction_type
        return reactions['SMARTS'][mask].values[0]
    except:
        raise ValueError(f'SMARTS of {reaction_type} not available.')


def get_reverse(
        const: str,
) -> str:
    return const[::-1].replace('{codon}'[::-1], '{codon}')


def get_complement(
        const: str,
) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in const)


def generate_const(
        const: str,
) -> str:
    seq = const['Sequence']
    if const['Reverse']:
        seq = get_reverse(seq)
    if const['Complement']:
        seq = get_complement(seq)
    return seq


def insert_codon(
        const: str,
        codon: str,
) -> str:
    return const.replace('{codon}', codon, 1)


def compute_product(
        smarts: str,
        smiles_1: str,
        smiles_2: str = None,
) -> str:
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    if smiles_2:
        react_1, react_2 = Chem.MolFromSmiles(smiles_1), Chem.MolFromSmiles(smiles_2)
        try:
            product = rxn.RunReactants((react_1, react_2))[0][0]
        except:
            product = rxn.RunReactants((react_2, react_1))[0][0]
    else:
        react = Chem.MolFromSmiles(smiles_1)
        product = rxn.RunReactant(react, 0)[0][0]
    return Chem.MolToSmiles(product)


def perform_reaction(
        reaction_type: str,
        reactions: pd.DataFrame,
        smiles_1: str,
        smiles_2: str,
) -> tp.Tuple:
    unimolecular = ['SR', 'DH', 'DT', 'SN2-1']
    if not reaction_type:
        return smiles_1
    smarts = get_smarts(reaction_type, reactions)
    if reaction_type in unimolecular:
        return compute_product(smarts, smiles_1)
    else:
        return compute_product(smarts, smiles_1, smiles_2)


def read_txt(
        path: str,
) -> tp.List:
    if str(path)[-2:] == 'gz':
        with gzip.open(path, 'rt') as file:
            return file.readlines()
    else:
        with open(path, 'r') as file:
            return file.readlines()


def write_txt(
        rows: tp.List,
        path: str = './',
        mode: str = 'a',
) -> None:
    with open(path, mode) as file:
        for row in rows:
            file.write('\t'.join(row))
            file.write('\n')


def read_json(
        path: str,
) -> tp.Dict:
    with open(path, 'r') as file:
        return json.load(file)


def write_json(
        data: tp.Dict,
        path: str,
) -> None:
    with open(path, 'w') as file:
        json.dump(data, file)




if __name__ == '__main__':

    ABF = '[CX3:1](=[O:2])[OX2;H1].[N;H2:4]>>[CX3:1](=[O:2])[N;H:4]'
    SR = '[#6:1][$([NX2-][NX2+]#[NX1]),$([NX2]=[NX2+]=[NX1-])]>>[#6:1][N;H2]'
    CuAAC = '[CX2:1]#[CX2;H1:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N-0:3][N-0:4]=[N-0:5]1'
    Suz = '[cX3:1][I].[#6:2][BX3]>>[cX3:1][#6:2]'
    Son = '[cX3:1][I].[CX2:2]#[CX2;H1:3]>>[cX3:1]-[CX2:3]#[CX2:2]'
    DH = '[cX3:1][I]>>[cX3;H1:1]'
    DT = '[#6:1][N;H2:2]>>[#6:1][N:2]=[N+]=[N-]'

    bb1 = 'C[C@H](N=[N+]=[N-])C(=O)N[C@@H](C)C(=O)NC(=O)[C@H](N)CS'
    bb2 = 'C#CCCC(O)=O'

    product = compute_product(CuAAC, bb1, bb2)
    print(product)

