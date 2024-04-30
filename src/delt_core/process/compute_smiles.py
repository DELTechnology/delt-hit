import os
from pathlib import Path
import typing as tp

import numpy as np

from .utils import get_smiles, generate_code, perform_reaction, write_txt, read_txt


def compute_smiles(
        libraries: tp.List,
        output_path: str,
) -> None:
    
    parent = output_path.parent

    # Compute SMILES.
    for i, library in enumerate(libraries, 1):
        perform_reaction_steps(i, library, parent)
    

    # Hybridize libraries.
    num_libs = len(libraries)

    for lib1 in range(1, num_libs + 1):
        for lib2 in range(1, num_libs + 1):
            if lib1 < lib2:
                hybridize(lib1, lib2, parent)
        
        if num_libs > 1:
            os.remove(parent / f'smiles_{lib1}.txt')


def perform_reaction_steps(
        index: int,
        library: tp.Tuple,
        output_path: str,
) -> None:
    
    tmp = output_path / 'tmp.txt'
    output_path = output_path / f'smiles_{index}.txt'

    header = [f'Scaffold_L{index}', f'BuildingBlock1_L{index}', f'Product_L{index}', 'Sequence']
    write_txt([header], output_path, 'w')

    bbs, scaffolds, reactions, consts = library

    # Perform reaction step 1.
    bb1 = bbs.pop(0)
    const1 = consts.iloc[0]
    rows = []

    for _, bb in bb1.iterrows():

        if bb['ReactionType']:  
            scaffold = get_smiles(bb['ScaffoldID'], scaffolds)
            product = scaffold
            
            for reaction_type in bb['ReactionType'].split(','):
                product = perform_reaction(reaction_type, reactions, bb['SMILES'], product)
        
        else:
            scaffold = 'None'
            product = bb['SMILES']
        
        code = generate_code(const1, bb['Codon'])
        rows += [[scaffold, bb['SMILES'], product, code]]
    
    write_txt(rows, output_path)
    
    # Perform reaction steps 2, ..., n.
    for i, bbn in enumerate(bbs, 2):
        header.insert(i, f'BuildingBlock{i}_L{index}')
        write_txt([header], tmp, 'w')
        idx = (consts['ID'] == i).idxmax()
        constn = consts.iloc[idx]
        
        with open(output_path, 'r') as interms:
            next(interms)

            for interm in interms:
                interm = interm.split()
                rows = []

                for _, bb in bbn.iterrows():
                    product = perform_reaction(bb['ReactionType'], reactions, interm[-2], bb['SMILES'])
                    code = interm[-1] + generate_code(constn, bb['Codon'])
                    rows += [[*interm[:-2], bb['SMILES'], product, code]]
                
                write_txt(rows, tmp)
                break
        
        os.rename(tmp, output_path)


def hybridize(
        lib1: int,
        lib2: int,
        output_path: str,
) -> None:
    
    output_file = output_path / f'smiles_{lib1}-{lib2}.txt'

    lib1 = output_path / f'smiles_{lib1}.txt'
    lib2 = output_path / f'smiles_{lib2}.txt'

    interms1 = read_txt(lib1)
    interms2 = read_txt(lib2)
    
    header1 = interms1.pop(0).split()[:-1]
    header2 = interms2.pop(0).split()[:-1]

    header = [*header1, *header2, "Sequence"]
    write_txt([header], output_file, 'w')
    
    for interm1 in interms1:
        interm1 = interm1.split()
        code1 = interm1.pop()
        
        rows = []

        for interm2 in interms2:
            interm2 = interm2.split()
            code = code1 + interm2.pop()

            rows += [[*interm1, *interm2, code]]
        
        write_txt(rows, output_file)

