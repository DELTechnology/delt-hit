import ast
import os
from pathlib import Path

import numpy as np

from .utils import load_data, get_smiles, generate_code, perform_reaction, write_to_txt


def compute_products(
    input_path: str,
    output_path: str = None,
) -> None:
    
    if not output_path:
        output_path = Path(input_path).parent / 'products.txt'
    tmp = Path(output_path).parent / 'tmp.txt'

    header = ['Scaffold', 'BuildingBlock1', 'Product', 'Sequence']
    write_to_txt(header, output_path, 'w')

    bbs, scaffolds, reactions, consts = load_data(input_path)

    # STEP 1

    bb1 = bbs.pop(0)
    const1 = consts.iloc[0]

    for i, bb in bb1.iterrows():
        scaffold = get_smiles(bb['ScaffoldID'], scaffolds)
        product = (scaffold,)
        
        for reaction_type in bb['ReactionType'].split(','):
            product = perform_reaction(reaction_type, reactions, bb['SMILES'], product[0])
        
        code = generate_code(const1, bb['Codon'])
        row = np.array([scaffold, bb['SMILES'], str(product), code])
        write_to_txt(row, output_path)
    
    # STEP N

    for i, bbn in enumerate(bbs, 2):
        header.insert(i, f'BuildingBlock{i}')
        write_to_txt(header, tmp, 'w')
        constn = consts[consts['ID'] ==  i]
        
        with open(output_path, 'r') as interms:
            next(interms)

            for interm in interms:
                interm = interm.split()
                interm_product = ast.literal_eval(interm[-2])[0]

                for _, bb in bbn.iterrows():
                    product = perform_reaction(bb['ReactionType'], reactions, interm_product, bb['SMILES'])
                    code = interm[-1] + generate_code(constn, bb['Codon'])
                    row = [*interm[:-2], bb['SMILES'], str(product), code]
                    write_to_txt(row, tmp)
        
        os.rename(tmp, output_path)

