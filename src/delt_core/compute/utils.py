import gzip
from typing import Optional

import numpy as np
import pandas as pd
from pydantic import BaseModel, validator, ValidationError
from rdkit import Chem
from rdkit.Chem import rdChemReactions


class BuildingBlock(BaseModel):
    ID: int
    SMILES: Optional[str] = None
    ScaffoldID: Optional[int] = None
    Codon: Optional[str] = None
    ReactionType: Optional[str] = None

    @validator('ReactionType')
    def set_reaction_type(cls, reaction_type):
        return reaction_type or ''


class Scaffold(BaseModel):
    ScaffoldID: Optional[int] = None
    SMILES: Optional[str] = None


class Reaction(BaseModel):
    ReactionType: str
    SMARTS: str


class ConstRegions(BaseModel):
    Sequence: str
    Reverse: bool
    Complement: bool


"""
class Library(BaseModel):
    BuildingBlocks: list[BuildingBlock]
    Scaffolds: dict[Scaffold]
    Reactions: dict[Reaction]
    ConstRegions: ConstRegions
"""


def validate_buliding_block(
        data: pd.DataFrame,
        model: BaseModel,
) -> list:
    bb = []
    records = data.to_dict(orient='records')
    for record in records:
        try:
            bb += [model(**record)]
        except ValidationError as e:
            print(f'Validation error: {e.json()}')
    return bb


def validate(
        data: pd.DataFrame,
        model: BaseModel,
) -> pd.DataFrame:
    records = data.to_dict(orient='records')
    for record in records:
        try:
            model(**record)
        except ValidationError as e:
            print(f'Validation error: {e.json()}')
    return data


def load_data(
        path: str,
) -> tuple[list, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scaffolds = validate(pd.read_excel(path, sheet_name='scaffolds'), Scaffold)
    reactions = validate(pd.read_excel(path, sheet_name='smarts'), Reaction)
    consts = validate(pd.read_excel(path, sheet_name='const'), ConstRegions)
    bbs = []
    try:
        while True:
            step = len(bbs) + 1
            bb = pd.read_excel(path, sheet_name=f'step{step}').replace({np.nan: None})
            bbs += [validate_buliding_block(bb, BuildingBlock)]
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
        const: pd.DataFrame,
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
) -> tuple:
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
) -> list:
    with open(path, 'r') as file:
        return file.readlines()


def write_txt(
        rows: list,
        path: str = './',
        mode: str = 'a',
) -> None:
    with open(path, mode) as file:
        for row in rows:
            file.write('\t'.join(row))
            file.write('\n')


def read_gzip(
        path: str,
) -> list:
    with gzip.open(path, 'rt') as file:
        return file.readlines()


def write_gzip(
        rows: list,
        path: str = './',
        mode: str = 'at',
) -> None:
    with gzip.open(path, mode) as file:
        for row in rows:
            file.write('\t'.join(row))
            file.write('\n')

