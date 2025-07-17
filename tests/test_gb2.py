# from delt_core.compute.utils import compute_product, perform_reaction
# from delt_core.compute.compute_smiles import perform_reaction, compute_smiles
from delt_core.cli.compute.cmds import compute_smiles
from pathlib import Path

library_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/GB2.xlsx')
compute_smiles(input_files=(library_path,))

# %%

from rdkit import Chem
from rdkit.Chem import AllChem

# Define the reaction
rxn = AllChem.ReactionFromSmarts('[CX3:1](=[O:2])[OX2;H1].[N;H2:4]>>[CX3:1](=[O:2])[N;H:4]')

# match secondary amin
rxn = AllChem.ReactionFromSmarts('[CX3:1](=[O:2])[OX2;H1].[N;H1:4]>>[CX3:1](=[O:2])[N:4]')

# Reactants
acid = Chem.MolFromSmiles('O=C(N)CCC(N=[N+]=[N-])C(O)=O')
amine = Chem.MolFromSmiles('C1CCNCC1')  # piperazine

# Run reaction
products = rxn.RunReactants((acid, amine))

# Show product SMILES
for i, product_tuple in enumerate(products):
    product = product_tuple[0]
    smiles = Chem.MolToSmiles(product)
    print(f"Product {i+1}: {smiles}")