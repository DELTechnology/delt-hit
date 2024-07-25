from collections import Counter
from collections.abc import Generator
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import QED

from .utils import write_gzip


PROPERTIES = ['ALERTS', 'ALOGP', 'AROM', 'HBA', 'HBD', 'MW', 'PSA', 'ROTB']
BINS = [None, 30, None, None, None, 30, 30, 20]


def read_gzip(
        file: Path,
        chunk_size: int = 1000,
) -> Generator:
    with gzip.open(file, 'rt') as f:
        while True:
            chunks = []
            try:
                for _ in range(chunk_size):
                    chunk = next(f)
                    chunks += [chunk]
            except StopIteration:
                if chunks:
                    yield chunks
                break
            yield chunks


def compute_qed_properties(
        mol_structures: list,
) -> list:
    return [QED.properties(mol) for mol in mol_structures]


def plot_properties(
        input_file: Path,
        output_dir: Path,
) -> None:

    for i, property in enumerate(PROPERTIES):
        
        plt.figure()
        output_file = output_dir / f'{property}.png'
        counter = Counter()

        chunks = read_gzip(input_file)
        for j, chunk in enumerate(chunks):
            if not j:
                chunk.pop(0)
            data = [float(line.split('\t')[i]) for line in chunk]
            counter.update(data)

        data = list(counter.elements())
        plt.hist(data, bins=BINS[i])
        plt.xlabel(property)
        plt.savefig(output_file, dpi=300)


def compute_properties(
        input_file: Path,
        index: int,
        output_file: Path = 'properties.txt.gz',
) -> None:

    chunks = read_gzip(input_file)
    for i, chunk in enumerate(chunks):
        
        if not i:
            chunk.pop(0)
            write_gzip([PROPERTIES], output_file, 'wt')
        
        smiles = [line.split('\t')[index] for line in chunk]
        mols = [Chem.MolFromSmiles(m) for m in smiles]
        qed_properties = compute_qed_properties(mols)
        properties = [[str(getattr(molecule, property)) for property in PROPERTIES] for molecule in qed_properties]
        write_gzip(properties, output_file)

