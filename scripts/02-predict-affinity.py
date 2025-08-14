import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# %%
counts_path = Path(
    '/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/evaluations/selection-4/143711a105ffe8e144c898e8114741df0cb28d5bb62dec6f22be4aa2846a1910.txt')
counts = pd.read_csv(counts_path, sep='\t')
# counts.sort_values('Count').tail(300).Code1.value_counts()
# code1 in {346, 246}
#         Count  Code1  Code2
# 98303    7102    246    126
# 138925   6924    346    126
# 139029   6418    346    230
# 139142   5776    346    344
# 98407    5401    246    230
counts = counts.set_index(['Code1', 'Code2'])

with open('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/NF2/index.pkl', 'rb') as f:
    index = pickle.load(f)

morgan_fps_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/NF2/morgan_fps.npz')
morgan_fps = sparse.load_npz(morgan_fps_path)
morgan_fps = morgan_fps.toarray()

bert_fps_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/NF2/bert_fps.npy')
bert_fps = np.load(bert_fps_path)

path = "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/NF_smiles.txt.gz"
df = pd.read_csv(path, compression='gzip', sep='\t')

df = df.rename(columns={'BuildingBlock1_L1': 'Code1', 'BuildingBlock2_L1': 'Code2'})
df = df[['Code1', 'Code2']]
# df = df.sample(1000)

morgan_fps = pd.DataFrame(morgan_fps)
morgan_fps = pd.concat([morgan_fps, df], axis=1, join='inner').set_index(['Code1', 'Code2'])
counts, morgan_fps = counts.align(morgan_fps, axis=0, join='inner')

bert_fps = pd.DataFrame(bert_fps)
bert_fps = pd.concat([bert_fps, df], axis=1, join='inner').set_index(['Code1', 'Code2'])
counts, bert_fps = counts.align(bert_fps, axis=0, join='inner')

# %%
import umap.plot

mapping = {0: (180 / 255, 188 / 255, 180 / 255), 346: (1, 0, 0), 246: (0, 1, 0)}

keep = np.random.rand((len(df))) > 0.8
subset_points = df.Code2.isin([346, 246]).values | keep

pdat = morgan_fps[subset_points]
labels = counts[keep].Code1.map(lambda x: {346: 346, 246: 246}.get(x, 0))
reducer = umap.UMAP(metric='jaccard').fit(X=pdat)

ax = umap.plot.points(reducer, background='black', labels=labels, color_key=mapping)
ax.figure.show()

# %%
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor
from sklearn.model_selection import cross_validate
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import make_scorer
from scipy.stats import spearmanr
import numpy as np


def spearman_corr(y_true, y_pred):
    return spearmanr(y_true, y_pred).correlation


spearman_scorer = make_scorer(spearman_corr, greater_is_better=True)

random_state = 0
scoring = {'mae': 'neg_mean_absolute_error', 'spearman': spearman_scorer},
models = {
    "RandomForest": RandomForestRegressor(
        n_estimators=50, max_depth=5, n_jobs=-1, random_state=random_state
    ),
    "XGBoost": XGBRegressor(
        n_estimators=50,
        max_depth=6,
        learning_rate=0.1,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="reg:squarederror",
        random_state=random_state,
    ),
    "kNN": KNeighborsRegressor(
        n_neighbors=5, weights="distance", p=2
    ),
    "MLP": MLPRegressor(
        hidden_layer_sizes=(128, 64),
        activation="relu",
        solver="adam",
        alpha=1e-4,
        batch_size="auto",
        learning_rate_init=1e-3,
        max_iter=200,
        random_state=random_state,
    ),
}

counts = counts.sample(1000)
y = counts.Count
X = counts[['Code1', 'Code2']]
for name, model in models.items():
    scores = cross_validate(estimator=model, X=X, y=y, scoring='neg_mean_absolute_error', cv=5, n_jobs=-1,
                            return_train_score=True)
    scores = pd.DataFrame(scores)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=random_state)

    model.fit(X_train, y_train)
    predictions = model.predict(X_test)
    corr, _ = spearmanr(predictions, y_test)
    mae = mean_absolute_error(y_test, predictions)
    print(f"{name:>12s} – Spearman correlation: {corr:.3f}, MAE: {mae:.5f}")

filter_ = (counts.Code1.isin([346, 246])) & (counts.Code2 == 1)
counts[filter_]
counts.groupby('Code1').Count.mean().sort_values()

counts = counts.set_index(['Code1', 'Code2']).squeeze()
y, X = counts.align(morgan_fps, join='inner', axis=0)
for name, model in models.items():
    scores = cross_validate(estimator=model, X=X, y=y, scoring='neg_mean_absolute_error', cv=5, n_jobs=-1,
                            return_train_score=True)
    scores = pd.DataFrame(scores)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=random_state)

    predictions = model.predict(X_test)
    corr, _ = spearmanr(predictions, y_test)
    mae = mean_absolute_error(y_test, predictions)
    print(f"{name:>14s} – Spearman correlation: {corr:.3f}, MAE: {mae:.5f}")

y, X = counts.align(bert_fps, join='inner', axis=0)
for name, model in models.items():
    scores = cross_validate(estimator=model, X=X, y=y, scoring=scoring, cv=5, n_jobs=-1,
                            return_train_score=True)
    scores = pd.DataFrame(scores)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=random_state)

    predictions = model.predict(X_test)
    corr, _ = spearmanr(predictions, y_test)
    mae = mean_absolute_error(y_test, predictions)
    print(f"{name:>14s} – Spearman correlation: {corr:.3f}, MAE: {mae:.5f}")

# %%
mean_estimator = counts.mean()
score = mean_absolute_error(counts, np.repeat(mean_estimator, len(counts)))
name = 'Mean Estimator'
print(f"{name:>124s} – Spearman correlation: {corr:.3f}, MAE: {mae:.5f}")
