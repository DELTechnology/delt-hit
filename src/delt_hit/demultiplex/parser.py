from pathlib import Path
import pandas as pd

def validate_config(config: dict):
    """Validate that the config is internally consistent.

    Args:
        config: Parsed configuration dictionary.

    Raises:
        AssertionError: If the library reactions are not all present in the catalog.
    """
    assert set(config['catalog']['reactions']) >= set(config['library']['reactions']), f'All reactions in library must be present in the `reactions` sheet'

def config_from_excel(path: Path):
    """Load a full configuration from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        The assembled configuration dictionary.
    """
    config = {}
    config['experiment'] = experiment_from_excel(path)
    config['selections'] = selections_from_excel(path)

    config['library'] = library_from_excel(path)
    config['catalog'] = catalog_from_excel(path)

    config['structure'] = structure_from_excel(path)
    # config['analyses'] = analyses_from_excel(path)
    config['whitelists'] = whitelists_from_excel(path)

    validate_config(config)
    return config

def library_from_excel(path: Path) -> dict:
    """Parse library topology metadata from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dictionary with library edges, nodes, and building block sheets.

    Raises:
        AssertionError: If required columns are missing values.
    """
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
    """Parse experiment metadata from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dict of experiment settings keyed by variable name.
    """
    experiment = pd.read_excel(path, sheet_name='experiment')
    return experiment.set_index('variable')['value'].to_dict()

def structure_from_excel(path: Path):
    """Parse structure metadata from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A list of structure entries as dictionaries.

    Raises:
        AssertionError: If structure names or types are invalid.
    """
    structure = pd.read_excel(path, sheet_name='structure')
    assert structure.name.str.match(r'^[SCB]').all(), "Structure `name` must start with 'S', 'B', or 'C' depending on type"
    assert structure.type.isin(['selection', 'building_block', 'constant']).all(), "Structure `type` must be one of 'selection', 'building_block', or 'constant'"
    return structure.to_dict('records')

def selections_from_excel(path: Path):
    """Parse selections metadata from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dict keyed by selection name with selection metadata.

    Raises:
        AssertionError: If selection names or primer combinations are invalid.
    """
    selections = pd.read_excel(path, sheet_name='selection')
    if 'date' in selections.columns:
        selections['date'] = pd.to_datetime(selections['date']).dt.strftime('%Y-%m-%d')
    selection_ids_to_name = get_selection_name_to_ids(path)
    assert selections.name.is_unique, "Selection `names` must be unique"
    selections['ids'] = list(map(list, selections['name'].map(selection_ids_to_name).tolist()))
    return selections.set_index('name').to_dict('index')

def analyses_from_excel(path: Path):
    """Group selections by analysis name from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dict mapping analysis name to selection names.
    """
    selections = pd.read_excel(path, sheet_name='selection')
    analyses = {}
    for grp, data in selections.groupby('analysis'):
        analyses[grp] = data.name.tolist()
    return analyses

def whitelists_from_excel(path: Path):
    """Parse codon whitelist sheets from an Excel workbook.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dict mapping whitelist names to codon records.

    Raises:
        AssertionError: If codons are missing or duplicated.
    """
    xf = pd.ExcelFile(path)
    sheets = set(xf.sheet_names)

    bbs_sheets = sorted(filter(lambda x: x.startswith('B'), sheets))

    selections = pd.read_excel(path, sheet_name='selection')
    selection_col_names = list(filter(lambda x: x.startswith('S'), selections.columns))

    constants = pd.read_excel(path, sheet_name='constant')
    assert constants.notna().any().any(), "`constant` cannot have empty cells"

    # %%
    whitelists = {}

    constants = constants.to_dict('records')
    for constant in constants:
        name, codon = constant.pop('name'), constant.pop('codon')
        whitelists[name] = [{'codon': codon}]

    for sheet in bbs_sheets:
        df = pd.read_excel(path, sheet_name=sheet)
        assert df.codon.nunique() == len(df), f"Codons for building blocks {sheet} must be unique"
        assert df.codon.notna().all(), f"Codons for building blocks {sheet} cannot be empty"

        df.index.name = 'index'
        df = df.reset_index()
        whitelists[sheet] = df.to_dict('records')

    assert len(set(selections[selection_col_names].to_records(index=False).tolist())) == len(selections), "Selection primers `S0` and `S1` must form unique combinations for each selection"
    for col_name in selection_col_names:
        df = selections[['name', col_name]].rename(columns={col_name: 'codon'})
        assert df.codon.notna().all(), f"Codons for selection primer {col_name} in `selection` cannot be empty"
        whitelists[col_name] = df.to_dict('records')

    return whitelists

def catalog_from_excel(path: Path):
    """Parse the catalog of compounds and reactions from Excel.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A dict containing compound and reaction catalogs.
    """
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
    """Build a selection name to primer ID mapping.

    Args:
        path: Path to the Excel configuration workbook.

    Returns:
        A mapping of selection names to primer ID tuples.

    Raises:
        AssertionError: If primer combinations are not unique.
    """
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
