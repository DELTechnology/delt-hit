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
    """
    Returns the SMILES of the respective scaffold.
    """
    try:
        mask = scaffolds['ScaffoldID'] == scaffold_id
        return scaffolds['SMILES'][mask].values[0]
    except:
        raise ValueError('SMILES not available.')


def get_smarts(
        reaction_type: str,
        reactions: pd.DataFrame,
) -> str:
    """
    Returns the SMARTS of the respective reaction.
    """
    try:
        mask = reactions['ReactionType'] == reaction_type
        return reactions['SMARTS'][mask].values[0]
    except:
        raise ValueError('SMARTS not available.')


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
    """
    Computes the product of the respective reaction.
    """
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    if smiles_2:
        reacts = (Chem.MolFromSmiles(smiles_1), Chem.MolFromSmiles(smiles_2))
        product = rxn.RunReactants(reacts)
    else:
        react = Chem.MolFromSmiles(smiles_1)
        product = rxn.RunReactant(react, 0)
    return Chem.MolToSmiles(product[0][0])


def perform_reaction(
        reaction_type: str,
        reactions: pd.DataFrame,
        smiles_1: str,
        smiles_2: str,
) -> tp.Tuple:
    """
    Performs the respective reaction.
    """
    if not reaction_type:
        return smiles_1
    smarts = get_smarts(reaction_type, reactions)
    if reaction_type == 'SR':
        return compute_product(smarts, smiles_2)
    elif reaction_type == 'DH':
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

    pass

