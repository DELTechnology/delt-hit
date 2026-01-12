from pathlib import Path
import pandas as pd

# path = Path('/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/template.xlsx')
# path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/251021_NF2_library_rechecked_enumerate_v1.xlsx')

def config_from_excel(path: Path):
    config = {}
    config['experiment'] = experiment_from_excel(path)
    config['selections'] = selections_from_excel(path)

    config['library'] = library_from_excel(path)
    config['catalog'] = catalog_from_excel(path)

    config['structure'] = structure_from_excel(path)
    # config['analyses'] = analyses_from_excel(path)
    config['whitelists'] = whitelists_from_excel(path)
    return config

def library_from_excel(path: Path) -> dict:
    rnx_g = pd.read_excel(path, sheet_name='reaction_graph')
    assert not rnx_g.educt_1.isna().any(), "All `educt_1` in `reaction_graph` must be filled"
    assert not rnx_g['product'].isna().any(), "All `product` in `reaction_graph` must be filled"
    assert not rnx_g.reaction.isna().any(), "All `reaction` in `reaction_graph` must be filled"

    educts = set()
    reactions = set()
    products = set()

    edges = [(i, j) for i, j in zip(rnx_g.educt_1, rnx_g.reaction)]
    edges += [(i, j) for i, j in zip(rnx_g.reaction, rnx_g['product'])]
    educts.update(rnx_g.educt_1)
    reactions.update(rnx_g.reaction)
    products.update(rnx_g['product'])

    rnx_g = rnx_g.dropna(subset=['educt_2'])
    edges += [(i, j) for i, j in zip(rnx_g.educt_2, rnx_g.reaction)]
    edges += [(i, j) for i, j in zip(rnx_g.reaction, rnx_g['product'])]
    educts.update(rnx_g.educt_2)

    xf = pd.ExcelFile(path)
    sheets = set(xf.sheet_names)
    bbs_sheets = sorted(filter(lambda x: x.startswith('B'), sheets))
    bb_edges = set()
    for sheet in bbs_sheets:
        df = pd.read_excel(path, sheet_name=sheet)
        df = df.astype(str)

        filter_ = df.smiles.notna()
        bb_edges.update([(sheet, r) for r in df.reaction[filter_]])
        bb_edges.update([(i, r) for i, r in zip(df.educt, df.reaction)])
        bb_edges.update([(r, p) for r,p in zip(df.reaction, df['product'])])

        products.update(df['product'])
        educts.update(df['educt'])
        reactions.update(df['reaction'])

    edges = [list(i) for i in edges]
    bb_edges = [list(i) for i in bb_edges]

    return dict(products=sorted(list(products)), educts=sorted(list(educts)), reactions=sorted(list(reactions)),
                bb_edges=sorted(list(bb_edges)), other_edges=sorted(list(edges)), building_blocks=sorted(list(bbs_sheets)))


def experiment_from_excel(path: Path):
    experiment = pd.read_excel(path, sheet_name='experiment')
    return experiment.set_index('variable')['value'].to_dict()

def structure_from_excel(path: Path):
    structure = pd.read_excel(path, sheet_name='structure')
    assert structure.name.str.match(r'^[SCB]').all(), "Structure `name` must start with 'S', 'B', or 'C' depending on type"
    assert structure.type.isin(['selection', 'building_block', 'constant']).all(), "Structure `type` must be one of 'selection', 'building_block', or 'constant'"
    return structure.to_dict('records')

def selections_from_excel(path: Path):
    selections = pd.read_excel(path, sheet_name='selection')
    if 'date' in selections.columns:
        selections['date'] = pd.to_datetime(selections['date']).dt.strftime('%Y-%m-%d')
    selection_ids_to_name = get_selection_name_to_ids(path)
    assert selections.name.is_unique, "Selection `names` must be unique"
    selections['ids'] = list(map(list, selections['name'].map(selection_ids_to_name).tolist()))
    return selections.set_index('name').to_dict('index')

def analyses_from_excel(path: Path):
    selections = pd.read_excel(path, sheet_name='selection')
    analyses = {}
    for grp, data in selections.groupby('analysis'):
        analyses[grp] = data.name.tolist()
    return analyses

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
        assert df.codon.nunique() == len(df), f"Codons for building blocks {sheet} must be unique"

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

    compounds = pd.read_excel(path, sheet_name='compounds')
    reactions = pd.read_excel(path, sheet_name='reactions')

    # %%
    catalog = {
        'compounds': compounds.set_index('name').to_dict('index'),
        'reactions': reactions.set_index('name').to_dict('index'),
    }

    return catalog

def get_selection_name_to_ids(path: Path)-> dict:
    df = pd.read_excel(path, sheet_name='selection')
    selection_col_names = list(filter(lambda x: x.startswith('S'), df.columns))

    assert len(df[selection_col_names].drop_duplicates()) == len(df), "S0, S1 combinations must be unique"

    df = df[['name', *selection_col_names]]

    for col_name in selection_col_names:
        mapping = df[[col_name]].drop_duplicates().reset_index().set_index(col_name)['index'].to_dict()
        df[col_name] = df[col_name].map(mapping)

    id_to_name = df.set_index(selection_col_names).name.to_dict()
    name_to_id = {v:k for k,v in id_to_name.items()}
    return name_to_id
