# from delt_core.compute.utils import compute_product, perform_reaction
# from delt_core.compute.compute_smiles import perform_reaction, compute_smiles
from delt_core.cli.compute.cmds import compute_smiles
from pathlib import Path
import pandas as pd

lib = pd.read_excel('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/GB2.xlsx', sheet_name='smarts')
library_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/GB2.xlsx')
output = compute_smiles(input_files=(library_path,))
# exit(0)

import pandas as pd
errors = pd.read_csv('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/GB2.err')
errors.error.value_counts()
no_product = errors[errors.error == 'NoProduct']
no_product.smiles_1.nunique()
no_product.smiles_2.nunique()
no_product.smarts.unique()

no_product = no_product.groupby('smarts').sample(25)
no_product.to_csv('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/no-products.csv')

# %%

rxn = rdChemReactions.ReactionFromSmarts(smarts)
react_1, react_2 = Chem.MolFromSmiles(smiles_1), Chem.MolFromSmiles(smiles_2)

rxn.RunReactants((react_1, react_2))
rxn.RunReactants((react_2, react_1))

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

# %%
import pandas as pd

df = pd.read_csv(
    "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/GB2_smiles.txt.gz",
    compression='gzip',
    sep='\t'  # or ',' or other, depending on actual delimiter
)

print(df.head())

pd.read_excel('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/GB2.xlsx', sheet_name='step1')
pd.read_excel('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/GB2.xlsx', sheet_name='step2')
