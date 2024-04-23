import ast
import os

import click
import numpy as np
import pandas as pd
import time

from .utils import load_data, get_smiles, generate_code, perform_reaction, write_to_txt

@click.command()
@click.argument('file', type=click.Path(writable=True, dir_okay=False))
@click.argument('num_steps', default=1, type=int)
def compute_products(
    file: str,
    num_steps: int = 1,
):
    cwd = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(cwd, '../../../data')
    input_path = os.path.join(data_path, file)
    output_path = os.path.join(data_path, 'products.txt')
    tmp = os.path.join(data_path, 'tmp.txt')

    header = ['Scaffold', 'BuildingBlock1', 'Product', 'Sequence']
    write_to_txt(header, output_path, 'w')

    bbs, scaffolds, reactions, consts = load_data(input_path, num_steps=num_steps)

    # STEP 1

    bb1 = bbs.pop(0)
    const1 = consts.iloc[0]

    for i, bb in bb1.iterrows():
        scaffold = get_smiles(bb['ScaffoldID'], scaffolds)
        product = scaffold
        
        for reaction_type in bb['ReactionType'].split(','):
            product = perform_reaction(reaction_type, reactions, bb['SMILES'], product)
        
        code = generate_code(const1, bb['Codon'])
        row = np.array([scaffold, bb['SMILES'], str(product), code])
        write_to_txt(row, output_path)
    
    # STEP N

    for i, bbn in enumerate(bbs, 2):
        header.insert(i, f'BuildingBlock{i}')
        write_to_txt(header, tmp, 'w')
        constn = consts.iloc[i-1]
        
        with open(output_path, 'r') as interms:
            next(interms)

            for _, interm in enumerate(interms):
                interm = interm.split()
                interm_product = ast.literal_eval(interm[-2])[0]

                for _, bb in bbn.iterrows():
                    product = perform_reaction(bb['ReactionType'], reactions, interm_product, bb['SMILES'])
                    code = interm[-1] + generate_code(constn, bb['Codon'])
                    row = [*interm[:-2], bb['SMILES'], str(product), code]
                    write_to_txt(row, tmp)
        
        os.rename(tmp, output_path)


if __name__ == '__main__':

    compute_products()

