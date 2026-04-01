# Chemistry Module Rules
- SMILES strings must be validated with RDKit before storage
- SMIRKS reaction templates must be tested with at least one example
- Building block SMILES should be canonical (RDKit.Chem.MolToSmiles)
- Never assume a reaction will succeed — always check for None products
