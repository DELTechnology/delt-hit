from pathlib import Path
import pandas as pd

# %%
counts_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/evaluations/selection-4/143711a105ffe8e144c898e8114741df0cb28d5bb62dec6f22be4aa2846a1910.txt')
counts = pd.read_csv(counts_path, sep='\t')

morgan_fp_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/NF2/morgan_fp.parquet')
morgan_fp = pd.read_parquet(morgan_fp_path, engine='fastparquet')

path = "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/NF_smiles.txt.gz"
df = pd.read_csv(path, compression='gzip', sep='\t')

df = df.rename(columns={'BuildingBlock1_L1': 'Code1', 'BuildingBlock2_L1': 'Code2'})
df = df[['Code1', 'Code2']]
morgan_fp = pd.concat([morgan_fp, df], axis=1).set_index(['Code1', 'Code2'])

# %%
counts = counts.set_index(['Code1', 'Code2']).squeeze()
counts, morgan_fp = counts.align(morgan_fp, join='inner', axis=0)

# %%
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_validate

reg = RandomForestRegressor(max_depth=5)
scores = cross_validate(estimator=reg, X=morgan_fp, y=counts, scoring='neg_mean_absolute_error', cv=5, n_jobs=-1)


