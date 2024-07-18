from typing import Optional

import pandas as pd
from pydantic import BaseModel, validator, ValidationError


class BuildingBlock(BaseModel):
    ID: int
    SMILES: Optional[str] = None
    ScaffoldID: Optional[int] = None
    Codon: Optional[str] = None
    ReactionType: Optional[str] = None

    @validator('ReactionType')
    def set_reaction_type(cls, reaction_type):
        return reaction_type or ''

    @validator('Codon')
    def set_codon(cls, codon):
        return codon or ''


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


