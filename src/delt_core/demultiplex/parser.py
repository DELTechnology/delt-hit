from pathlib import Path

import pandas as pd


def config_from_excel(path: Path):
    config = {}
    config['experiment'] = experiment_from_excel(path)
    config['structure'] = structure_from_excel(path)
    config['selections'] = selections_from_excel(path)
    config['catalog'] = catalog_from_excel(path)
    config['whitelists'] = whitelists_from_excel(path)
    return config

def experiment_from_excel(path: Path):
    experiment = pd.read_excel(path, sheet_name='experiment')
    return experiment.set_index('variable')['value'].to_dict()

def structure_from_excel(path: Path):
    structure = pd.read_excel(path, sheet_name='structure')
    return structure.to_dict('records')

def selections_from_excel(path: Path):
    selections = pd.read_excel(path, sheet_name='selection')
    selections['date'] = pd.to_datetime(selections['date']).dt.strftime('%Y-%m-%d')
    selection_ids_to_name = get_selection_name_to_ids(path)
    assert selections.name.is_unique
    selections['ids'] = selections['name'].map(selection_ids_to_name)
    return selections.set_index('name').to_dict('index')

def whitelists_from_excel(path: Path):
    xf = pd.ExcelFile(path)
    sheets = set(xf.sheet_names)

    bbs_sheets = sorted(filter(lambda x: x.startswith('B'), sheets))

    selections = pd.read_excel(path, sheet_name='selection')
    selection_col_names = list(filter(lambda x: x.startswith('S'), selections.columns))

    constants = pd.read_excel(path, sheet_name='constant')

    # %%
    whitelists = {}

    constants = constants.to_dict('records')
    for constant in constants:
        name, codon = constant.pop('name'), constant.pop('codon')
        whitelists[name] = [{'codon': codon}]

    for sheet in bbs_sheets:
        df = pd.read_excel(path, sheet_name=sheet)
        df.index.name = 'index'
        df = df.reset_index()
        whitelists[sheet] = df.to_dict('records')

    for col_name in selection_col_names:
        df = selections[['name', col_name]].rename(columns={col_name: 'codon'})
        whitelists[col_name] = df.to_dict('records')

    return whitelists


def catalog_from_excel(path: Path):
    xf = pd.ExcelFile(path)
    sheets = set(xf.sheet_names)

    scaffolds = pd.read_excel(path, sheet_name='scaffolds')
    reactions = pd.read_excel(path, sheet_name='reactions')

    # %%
    catalog = {
        'scaffolds': scaffolds.set_index('name').to_dict('index'),
        'reactions': reactions.set_index('name').to_dict('index'),
    }

    return catalog

def get_selection_name_to_ids(path: Path)-> dict:
    df = pd.read_excel(path, sheet_name='selection')
    selection_col_names = list(filter(lambda x: x.startswith('S'), df.columns))

    df = df[['name', *selection_col_names]]

    for col_name in selection_col_names:
        mapping = df[[col_name]].drop_duplicates().reset_index().set_index(col_name)['index'].to_dict()
        df[col_name] = df[col_name].map(mapping)

    id_to_name = df.set_index(selection_col_names).name.to_dict()
    name_to_id = {v:k for k,v in id_to_name.items()}
    return name_to_id
