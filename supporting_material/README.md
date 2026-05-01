# Supporting Material

## Example Workflow

The files below are the data used in the paper to demonstrate the full set of DELT-Hit capabilities.

Download the files:

- `template.xlsx`
- `campaign.fastq.gz`
- `analysis.yaml`

from [figshare](https://doi.org/10.6084/m9.figshare.31198468).

## Publication Reproductions

This directory contains rerunnable reproductions for the following publications:

- Favalli et al. DOI: [10.1038/s41557-021-00660-y](https://doi.org/10.1038/s41557-021-00660-y)
- Pure-DEL DOI: [10.1126/science.adn3412](https://doi.org/10.1126/science.adn3412)

The rerunnable scripts are:

- [supporting_material/experiments/favalli/run.sh](experiments/favalli/run.sh)
- [supporting_material/experiments/pure-del/run.sh](experiments/pure-del/run.sh)

Please adjust the `cd` command in each script to match your local folder structure before running.

Exact-match validation of the published and DELT-Hit selection counts can be performed with:

- [supporting_material/experiments/favalli/compare_selection.py](experiments/favalli/compare_selection.py)
- [supporting_material/experiments/pure-del/compare_selection.py](experiments/pure-del/compare_selection.py)

Example usage:

```bash
uv run python supporting_material/experiments/favalli/compare_selection.py
uv run python supporting_material/experiments/pure-del/compare_selection.py
```

## Open Questions

### Pure-DEL
- Counts for Lane 1/2 are identical, why?