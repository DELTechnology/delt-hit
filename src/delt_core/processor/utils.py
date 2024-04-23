import typing as tp

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions


def load_data(
    path: str,
    num_steps: int = 1,
):
    """
    Returns the building blocks, the scaffolds, and the reactions.
    """
    scaffolds = pd.read_excel(path, sheet_name='scaffolds')
    reactions = pd.read_excel(path, sheet_name='smarts')
    consts = pd.read_excel(path, sheet_name='const')
    bbs = []
    try:
        for step in range(1, num_steps+1):
            bbs += [pd.read_excel(path, sheet_name=f'step{step}')]
        return bbs, scaffolds, reactions, consts
    except:
        raise ValueError('Adjust the number of reaction steps.')


def get_smiles(
    scaffold_id: str,
    scaffolds: pd.DataFrame,
):
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
):
    """
    Returns the SMARTS of the respective reaction.
    """
    try:
        mask = reactions['ReactionType'] == reaction_type
        return reactions['SMARTS'][mask].values[0]
    except:
        raise ValueError('SMARTS not available.')


def get_reverse(
    codon: str,
):
    return codon[::-1]


def get_complement(
    codon: str,
):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    codon = ''.join(complement.get(base, base) for base in codon)
    return codon


def generate_code(
    const: str,
    codon: str,
):
    if const['Reverse']:
        codon = get_reverse(codon)
    if const['Complement']:
        codon = get_complement(codon)
    return const['Sequence'].replace('{codon}', codon)


def compute_product(
    smarts: str,
    smiles_1: str,
    smiles_2: str = None,
):
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
):
    """
    Performs the respective reaction.
    """
    if reaction_type == '-':
        return (smiles_1,)
    if reaction_type == 'H':
        return (smiles_1, smiles_2)
    smarts = get_smarts(reaction_type, reactions)
    if reaction_type == 'SR':
        return compute_product(smarts, smiles_2)
    elif reaction_type == 'DH':
        return (compute_product(smarts, smiles_1),)
    else:
        return (compute_product(smarts, smiles_1, smiles_2),)


def write_to_txt(
    product: tp.List,
    path: str = './',
    mode: str = 'a',
):
    with open(path, mode) as file:
        file.write('\t'.join(product))
        # file.write('\t'.join(map(str, product)))
        file.write('\n')

