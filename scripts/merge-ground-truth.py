import pandas as pd

nf = pd.read_excel('/Users/adrianomartinelli/projects/delt/delt-core/paper/NF.xlsx', sheet_name='B0')
truth = pd.read_excel('/Users/adrianomartinelli/projects/delt/delt-core/paper/20231011_NF_library_Step1and2_withSmiles_CheckedAndCorrect.xlsx', sheet_name='step1')
assert all(nf.codon == truth.codons)

# nf.loc[nf.smiles.isna(), 'Structure'] = ['NaN1', 'NaN2', 'NaN3']
smiles_to_reaction = nf.set_index(['smiles', 'reactant'])['reaction'].to_dict()
smiles_to_reactant = nf.set_index('smiles')['reactant'].to_dict()

b0 = pd.merge(nf, truth, left_on='codon', right_on='codons', how='outer', suffixes=('_nf', '_truth'))
# b0 = b0.rename(columns=dict(smiles='smiles_nf'))
# b0['smiles'] = b0.Structure
# b0['reaction'] = b0.smiles.map(smiles_to_reaction)
# b0['reactant'] = b0.smiles.map(smiles_to_reactant)
# b0 = b0[['CdId', 'smiles', 'codon', 'reaction', 'reactant', 'product']]

nf = pd.read_excel('/Users/adrianomartinelli/projects/delt/delt-core/paper/NF.xlsx', sheet_name='B1')
truth = pd.read_excel('/Users/adrianomartinelli/projects/delt/delt-core/paper/20231011_NF_library_Step1and2_withSmiles_CheckedAndCorrect.xlsx', sheet_name='step 2')
assert all(nf.codon == truth[:-1].codon)

# smiles_to_reaction = nf.set_index('smiles')['reaction'].to_dict()
# smiles_to_reactant = nf.set_index('smiles')['reactant'].to_dict()

b1 = pd.merge(nf, truth, left_on='codon', right_on='codon', how='outer', suffixes=('_nf', '_truth'))
# b1 = b1.rename(columns=dict(smiles_truth='smiles'))
# b1['reaction'] = b0.smiles.map(smiles_to_reaction)
# b1['reactant'] = b0.smiles.map(smiles_to_reactant)
# b1 = b1[['CdId', 'smiles', 'codon', 'reaction', 'reactant', 'product']]

b0 = b0.sort_values(by=['CdId'])
# b0.loc[b0.Structure == 'N', 'smiles'] = pd.NA

# b1.loc[b1.smiles_truth == '[H]', 'smiles_truth'] = pd.NA
# b1 = b1[b1.smiles_truth != '[I]']
b1 = b1.sort_values(by=['CdId'])

out_path = '/Users/adrianomartinelli/projects/delt/delt-core/paper/NF_merged.xlsx'
with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
    b0.to_excel(writer, sheet_name='B0', index=False)
    b1.to_excel(writer, sheet_name='B1', index=False)