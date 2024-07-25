from pathlib import Path


def init(
        root: Path = None,
) -> None:
    if not root:
        root = Path.cwd()
    dirs = ['libraries', 'fastq_files', 'selections']
    for dir in dirs:
        path = Path(root) / dir
        path.mkdir(exist_ok=True)

