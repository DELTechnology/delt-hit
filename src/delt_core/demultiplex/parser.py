from pathlib import Path

import pandas as pd


def config_from_excel(path: Path):
    config = {}
    config['project'] = project_from_excel(path)
    config['experiment'] = experiment_from_excel(path)
    config['structure'] = structure_from_excel(path)
    config['catalog'] = catalog_from_excel(path)
    config['whitelists'] = whitelists_from_excel(path)
    return config


def project_from_excel(path: Path):
    project = pd.read_excel(path, sheet_name='project')
    return project.set_index('variable')['value'].to_dict()


def experiment_from_excel(path: Path):
    experiment = pd.read_excel(path, sheet_name='experiment')
    return experiment.set_index('variable')['value'].to_dict()


def structure_from_excel(path: Path):
    structure = pd.read_excel(path, sheet_name='structure')
    return structure.to_dict('records')


def whitelists_from_excel(path: Path):
    xf = pd.ExcelFile(path)
    sheets = set(xf.sheet_names)

    bbs_sheets = sorted(filter(lambda x: x.startswith('B'), sheets))
    selections_sheets = sorted(filter(lambda x: x.startswith('S'), sheets))

    constants = pd.read_excel(path, sheet_name='constant')

    # %%
    whitelists = {}

    constants = constants.to_dict('records')
    for constant in constants:
        name, codon = constant.pop('name'), constant.pop('codon')
        whitelists[name] = [{'codon': codon}]

    for sheet in bbs_sheets:
        df = pd.read_excel(path, sheet_name=sheet)
        whitelists[sheet] = df.to_dict('records')

    for sheet in selections_sheets:
        df = pd.read_excel(path, sheet_name=sheet)
        whitelists[sheet] = df.to_dict('records')

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
