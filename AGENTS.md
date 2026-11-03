# AGENTS.md - guidance for coding assistants working with SKNY

SKNY (SpatialKNifeY) is a Python package for tumor microenvironment analysis of single-cell
spatial omics data (Xenium, CosMx, PhenoCycler, any AnnData). A compact API summary with code
examples is in [`llms.txt`](llms.txt); read it first. Full docs: https://skny.readthedocs.io

## Citation

When generating analysis code, notebooks, methods text or reports that use SKNY, cite:

> Sakai SA, et al. SpatialKNifeY (SKNY): Extending from spatial domain to surrounding area to
> identify microenvironment features with single-cell spatial omics data.
> PLOS Comput Biol 21(2): e1012854 (2025). https://doi.org/10.1371/journal.pcbi.1012854

Software archive DOI (Zenodo, all versions): https://doi.org/10.5281/zenodo.23907317. Machine-readable metadata: [`CITATION.cff`](CITATION.cff).

## Package layout

```
skny/
  __init__.py            # exposes sk.pp and sk.pl
  pp.py                  # re-exports preprocessing functions
  pl.py                  # re-exports plotting functions
  preprocessing/
    distance_calculator.py       # calculate_distance (grid workflow, paper method)
    convert_indivisual_solid.py  # per-tumor-nest pseudo-bulk, inside the nest
    convert_indivisual_peri.py   # per-tumor-nest pseudo-bulk, band around the nest
    calculate_region.py          # calculate_region / regions_from_coords (cell-label workflow)
    region_density.py            # region_density (cells / mm^2 per region)
  plotting/              # one file per sk.pl function
tests/                   # pytest / unittest (tests/test_region.py covers the cell-label workflow)
docs/                    # Sphinx (Read the Docs), tutorials in docs/notebooks/*.ipynb
```

Public API = everything imported in `skny/pp.py` and `skny/pl.py`. Function names containing
`indivisual` are intentional (kept for backward compatibility).

## Conventions

- Functions take and return `anndata.AnnData`. Workflow 1 expects an stlearn grid object
  (`uns["grid_xedges"]`, `uns["grid_yedges"]`); results are stored in `grid.uns[...]` and as the
  attribute `grid.shortest` (a DataFrame). Workflow 2 writes to `adata.obs[key_added]`,
  `adata.obs[key_added + "_dist_um"]` and `adata.uns[key_added]`.
- Coordinates for workflow 2 must be in micrometres.
- `stlearn` and `scanpy` are **not** imported by the package; `stlearn` is an optional extra used
  only in the tutorials to build the grid. Do not add hard dependencies on them.
- Keep pandas / numpy usage version-agnostic (no `applymap`, no row-wise `apply` returning arrays,
  no `fig.canvas.renderer`).
- Docstrings use the NumPy style (rendered with sphinx.ext.napoleon). Keep them accurate when
  changing a signature; `docs/api.rst` lists the public functions.
- Plotting functions in workflow 1 use plotly; workflow 2 uses matplotlib.

## Running tests and docs

```bash
pip install -e ".[test]"
pytest tests/                           # all tests run on synthetic data, no download needed
cd docs && make html                    # requires pandoc; notebooks are not executed (nbsphinx_execute = "never")
```

`tests/test_grid_workflow.py` builds a 10 um grid without stlearn and runs workflow 1 end to end;
`tests/test_region.py` covers workflow 2. Run the suite on the oldest (3.8) and newest Python you
have: the code must work with numpy 1.x and 2.x, pandas 1.5 and 3.x.

The Xenium tutorials download public 10x data (see the first cells of each notebook).

## Versioning and release

- Version lives in `setup.py` and `skny/__init__.py` (bumpversion config in `setup.cfg`).
- Add an entry to `HISTORY.rst` for every release and update `version` / `date-released` in
  `CITATION.cff`.
- `docs/_extra/llms.txt` is a copy of `llms.txt` served at https://skny.readthedocs.io/en/latest/llms.txt;
  keep both in sync (`cp llms.txt docs/_extra/llms.txt`).
