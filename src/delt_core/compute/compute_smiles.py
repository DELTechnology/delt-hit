import gzip
import os
from pathlib import Path
import tempfile

import pandas as pd

from .utils import (
    get_smiles,
    generate_const,
    insert_codon,
    perform_reaction,
    read_gzip,
    write_gzip,
)


def compute_smiles(
        libraries: list,
        output_path: str,
) -> None:
    
    parent = output_path.parent
    with tempfile.TemporaryDirectory(dir=parent) as tmp:
        tmp = Path(tmp)

        for i, library in enumerate(libraries, 1):  
            perform_reaction_steps(i, library, tmp)
        
        if len(libraries) > 1:
            hybridize(tmp)
            merge_excel_files(libraries, 'hybrid.xlsx')
    
        os.rename(tmp / 'smiles1.txt', output_path)


def perform_reaction_steps(
        index: int,
        library: tuple,
        output_path: str,
        write_indices: bool = True,
) -> None:

    tmp = output_path / 'tmp.txt'
    output_path = output_path / f'smiles{index}.txt'

    header = [f'Scaffold_L{index}', f'BuildingBlock1_L{index}', f'Product_L{index}', f'Sequence_L{index}']
    write_gzip([header], output_path, 'wt')

    bbs, scaffolds, reactions, consts = library
    const = generate_const(consts.iloc[0])

    # Perform reaction step 1.
    bb1, *bbs = bbs
    rows = []

    for _, bb in bb1.iterrows():

        if bb['ScaffoldID']:  
            scaffold_id = bb['ScaffoldID']
            scaffold = get_smiles(scaffold_id, scaffolds)
            product = scaffold
            
            if bb['SMILES']:
                for reaction_type in bb['ReactionType'].split(','):
                    product = perform_reaction(reaction_type.strip(), reactions, product, bb['SMILES'])
        
        else:
            scaffold = scaffold_id = 'None'
            product = bb['SMILES']
        
        code = insert_codon(const, bb['Codon'])
        if write_indices:
            rows += [[str(scaffold_id), str(bb['ID']), product, code]]
        else:
            rows += [[scaffold, bb['SMILES'], product, code]]
        break
    
    write_gzip(rows, output_path)
    
    # Perform reaction steps 2, ..., n.
    for i, bbn in enumerate(bbs, 2):
        header.insert(i, f'BuildingBlock{i}_L{index}')
        write_gzip([header], tmp, 'wt')
        
        with gzip.open(output_path, 'rt') as interms:
            next(interms)

            for interm in interms:
                interm = interm.split()
                rows = []

                for _, bb in bbn.iterrows():
                    product = interm[i]
                    if bb['SMILES']:
                        for reaction_type in bb['ReactionType'].split(','):
                            product = perform_reaction(reaction_type.strip(), reactions, product, bb['SMILES'])
                    
                    code = insert_codon(interm[-1], bb['Codon'])
                    if write_indices:
                        rows += [[*interm[:i], str(bb['ID']), product, code]]
                    else:
                        rows += [[*interm[:i], bb['SMILES'], product, code]]
                
                write_gzip(rows, tmp)
        
        os.rename(tmp, output_path)


def hybridize(
        output_path: str,
) -> None:
    
    output_file = output_path / 'smiles1.txt'
    tmp = output_path / 'tmp.txt'

    interms1 = read_gzip(output_file)
    header1 = interms1.pop(0).split()[:-1]

    for i, interm1 in enumerate(interms1):
        interm1 = interm1.split()
        code1 = interm1.pop()
        rows = []
        
        path = output_path / f'smiles2.txt'
        interms2 = read_gzip(path)
        header2 = interms2.pop(0).split()[:-1]
        
        if not i:
            header = [*header1, *header2, "Sequence"]
            write_gzip([header], tmp, 'wt')

        for interm2 in interms2:
            interm2 = interm2.split()
            code = code1 + interm2.pop()
            rows += [[*interm1, *interm2, code]]
    
        write_gzip(rows, tmp)
    
    os.rename(tmp, output_file)
    os.remove(path)


def merge_excel_files(
        libraries: list,
        output_file: Path,
) -> None:
    l1, l2 = libraries
    names = ['step', 'scaffolds', 'smarts', 'const']
    with pd.ExcelWriter(output_file) as writer:
        for sheet_l1, sheet_l2, name in zip(l1, l2, names):
            if name == 'step':
                step = 0
                for bbs in [sheet_l1, sheet_l2]:
                    for bb in bbs:
                        step += 1
                        bb.to_excel(writer, sheet_name=f'{name}{step}', index=False)
            elif name == 'const':
                sequence = generate_const(sheet_l1.iloc[0]) + generate_const(sheet_l2.iloc[0])
                const = {
                    'Sequence': [sequence],
                    'Reverse': [0],
                    'Complement': [0],
                }
                pd.DataFrame(const).to_excel(writer, sheet_name=name, index=False)
            else:
                pd.DataFrame().to_excel(writer, sheet_name=name, index=False)

