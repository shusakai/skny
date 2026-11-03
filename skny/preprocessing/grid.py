"""Gridding of single-cell spatial data and a light-weight Xenium reader.

:func:`grid` reproduces the output format of ``stlearn.tl.cci.grid`` (bins named ``grid_<n>``
with ``n = col * n_row + row`` starting at 0, ``uns["grid_xedges"]`` / ``uns["grid_yedges"]``)
without the stlearn / numba dependency, so that the expression-based workflow of SKNY
(:func:`skny.pp.calculate_distance`) can be run with ``pip install skny`` alone.
"""
import os
import glob
import warnings

import numpy as np
import pandas as pd
import anndata as ad


# --------------------------------------------------------------------------------------------
# helpers shared by the grid workflow
# --------------------------------------------------------------------------------------------
def grid_shape(grid):
    """Return ``(n_row, n_col)`` of a gridded AnnData (from ``uns["grid_yedges"]`` / ``uns["grid_xedges"]``)."""
    return len(grid.uns["grid_yedges"]) - 1, len(grid.uns["grid_xedges"]) - 1


def bin_names(grid):
    """Names of *all* bins of a gridded AnnData in column-major order (``n = col * n_row + row``).

    ``stlearn.tl.cci.grid`` and :func:`skny.pp.grid` number the bins from 0 (``grid_0`` ...);
    SKNY <= 0.3.1 assumed numbering from 1, which shifted every mask by one bin. The convention is
    detected (``uns["grid_bin_offset"]`` written by :func:`skny.pp.grid`, otherwise from the bin
    centre coordinates in ``obs``, otherwise from the smallest / largest bin number; default 0)
    so that objects of both kinds are handled correctly.
    """
    n_row, n_col = grid_shape(grid)
    return [f"grid_{i + _bin_offset(grid)}" for i in range(n_row * n_col)]


def _bin_offset(grid):
    if "grid_bin_offset" in grid.uns:
        return int(grid.uns["grid_bin_offset"])
    names = [n for n in grid.obs_names if str(n).startswith("grid_")]
    if not names:
        return 0
    nums = np.array([int(str(n).split("_", 1)[1]) for n in names], dtype=int)
    n_row, n_col = grid_shape(grid)
    if "imagecol" in grid.obs and "imagerow" in grid.obs:
        xe, ye = np.asarray(grid.uns["grid_xedges"]), np.asarray(grid.uns["grid_yedges"])
        k = min(len(names), 50)
        xs = grid.obs.loc[names[:k], "imagecol"].to_numpy(dtype=float)
        ys = grid.obs.loc[names[:k], "imagerow"].to_numpy(dtype=float)
        col = np.clip(np.searchsorted(xe, xs, side="right") - 1, 0, n_col - 1)
        row = np.clip(np.searchsorted(ye, ys, side="right") - 1, 0, n_row - 1)
        off = nums[:k] - (col * n_row + row)
        if (off == off[0]).all() and off[0] in (0, 1):
            return int(off[0])
    if nums.min() == 0:
        return 0
    if nums.max() == n_row * n_col:
        return 1
    return 0


# --------------------------------------------------------------------------------------------
# gridding
# --------------------------------------------------------------------------------------------
def grid(
    adata, bin_um=10.0, x_key="imagecol", y_key="imagerow", spatial_key=None,
    n_row=None, n_col=None, use_label=None, verbose=True,
):
    """Sum the expression of all cells falling into square bins (default 10 um).

    Drop-in replacement for ``stlearn.tl.cci.grid`` (same bin numbering, same ``uns`` keys)
    that only needs numpy / pandas / anndata. The output is the input of
    :func:`skny.pp.calculate_distance`.

    Parameters
    ----------
    adata : AnnData
        Cell-level object. Coordinates are read from ``obs[x_key]`` / ``obs[y_key]`` (as set by
        :func:`skny.pp.read_xenium` or stlearn readers) or, if ``spatial_key`` is given, from
        ``obsm[spatial_key][:, :2]``. The unit of the coordinates is the unit of ``bin_um``
        (micrometres for Xenium / CosMx / PhenoCycler).
    bin_um : float
        Side length of a bin. Ignored when both ``n_row`` and ``n_col`` are given.
    n_row, n_col : int or None
        Number of rows / columns of the grid (as in stlearn). Default: derived from ``bin_um``.
    use_label : str or None
        ``obs`` column with cell labels. If given, the fraction of each label per bin is stored
        in ``uns[use_label]`` and the majority label in ``obs[use_label]``.
    verbose : bool

    Returns
    -------
    AnnData
        One observation per non-empty bin, named ``grid_<n>`` with ``n = col * n_row + row``
        (0-based). ``obs["imagecol"]`` / ``obs["imagerow"]`` = bin centre, ``obs["n_cells"]``,
        ``obsm["spatial"]``, ``uns["grid_xedges"]``, ``uns["grid_yedges"]``, ``uns["grid_counts"]``
        (n_col x n_row matrix of cell counts). ``uns["spatial"]`` is copied from ``adata`` if present.

    Examples
    --------
    >>> adata = sk.pp.read_xenium("xenium_output_dir/")
    >>> grid = sk.pp.grid(adata, bin_um=10)
    >>> grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"])
    """
    if spatial_key is not None:
        xy = np.asarray(adata.obsm[spatial_key], dtype=float)[:, :2]
        xs, ys = xy[:, 0], xy[:, 1]
    else:
        if x_key not in adata.obs or y_key not in adata.obs:
            raise KeyError(
                f"obs['{x_key}'] / obs['{y_key}'] not found; pass spatial_key='spatial' to use obsm, "
                "or load the data with skny.pp.read_xenium"
            )
        xs = adata.obs[x_key].to_numpy(dtype=float)
        ys = adata.obs[y_key].to_numpy(dtype=float)
    if not (np.isfinite(xs).all() and np.isfinite(ys).all()):
        raise ValueError("coordinates contain NaN / inf")

    if n_row is None or n_col is None:
        if n_col is None:
            n_col = max(int((xs.max() - xs.min()) / bin_um), 1)
        if n_row is None:
            n_row = max(int((ys.max() - ys.min()) / bin_um), 1)
    n_row, n_col = int(n_row), int(n_col)
    if verbose:
        print(f"Gridding {adata.n_obs} cells into {n_row} rows x {n_col} cols ...")

    xedges = np.linspace(xs.min(), xs.max(), n_col + 1)
    yedges = np.linspace(ys.min(), ys.max(), n_row + 1)
    col = np.clip(np.searchsorted(xedges, xs, side="right") - 1, 0, n_col - 1)
    row = np.clip(np.searchsorted(yedges, ys, side="right") - 1, 0, n_row - 1)
    n = col * n_row + row  # stlearn numbering
    n_squares = n_row * n_col

    # expression sums per bin (sparse-safe)
    X = adata.X
    import scipy.sparse as sp
    M = sp.csr_matrix((np.ones(adata.n_obs), (n, np.arange(adata.n_obs))), shape=(n_squares, adata.n_obs))
    grid_expr = M @ (X if sp.issparse(X) else np.asarray(X))
    if sp.issparse(grid_expr):
        grid_expr = grid_expr.toarray()
    grid_expr = np.asarray(grid_expr, dtype=np.float32)

    counts = np.bincount(n, minlength=n_squares)
    cx = (xedges[:-1] + xedges[1:]) / 2
    cy = (yedges[:-1] + yedges[1:]) / 2
    ii = np.arange(n_squares)
    coords = np.c_[cx[ii // n_row], cy[ii % n_row]]

    g = ad.AnnData(
        grid_expr,
        obs=pd.DataFrame({"imagecol": coords[:, 0], "imagerow": coords[:, 1], "n_cells": counts},
                         index=[f"grid_{i}" for i in ii]),
        var=adata.var.copy(),
    )
    g.obsm["spatial"] = coords
    if "spatial" in adata.uns:
        g.uns["spatial"] = adata.uns["spatial"]

    if use_label is not None:
        labels = adata.obs[use_label].astype(str).to_numpy()
        cats = np.unique(labels)
        L = sp.csr_matrix((np.ones(adata.n_obs), (n, np.searchsorted(cats, labels))), shape=(n_squares, len(cats))).toarray()
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = L / L.sum(1, keepdims=True)
        info = pd.DataFrame(frac, index=g.obs_names, columns=cats)
        g.obs[use_label] = pd.Categorical(np.where(counts > 0, cats[np.argmax(L, 1)], "nan"))
        g.uns[use_label] = info

    keep = counts > 0
    g = g[keep].copy()
    if use_label is not None:
        g.uns[use_label] = g.uns[use_label].loc[g.obs_names]
    g.uns["grid_counts"] = counts.reshape(n_col, n_row)
    g.uns["grid_bin_offset"] = 0
    g.uns["grid_xedges"] = xedges.astype(float)
    g.uns["grid_yedges"] = yedges.astype(float)
    return g


# --------------------------------------------------------------------------------------------
# Xenium reader
# --------------------------------------------------------------------------------------------
def _read_10x_h5(path):
    """Minimal 10x Genomics HDF5 feature-barcode matrix reader (cells x features, CSR)."""
    import h5py
    import scipy.sparse as sp

    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = tuple(int(i) for i in m["shape"][:])  # (n_features, n_barcodes)
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape).T.tocsr()
        barcodes = m["barcodes"][:].astype(str)
        feat = m["features"]
        var = pd.DataFrame(index=feat["name"][:].astype(str))
        var["gene_ids"] = feat["id"][:].astype(str)
        if "feature_type" in feat:
            var["feature_types"] = feat["feature_type"][:].astype(str)
    return X, barcodes, var


def read_xenium(path, cells_file=None, library_id=None, genes_only=True, verbose=True):
    """Read a 10x Xenium run (cell x gene counts + cell centroids) into an AnnData for SKNY.

    Replacement for ``stlearn.ReadXenium`` that needs only ``h5py`` (installed with anndata).

    Parameters
    ----------
    path : str
        Either the Xenium output directory (containing ``cell_feature_matrix.h5`` and
        ``cells.csv.gz`` / ``cells.parquet``) or the path of the ``*cell_feature_matrix.h5`` file
        (then ``cells_file`` must be given).
    cells_file : str or None
        Path of the cell summary table (``*cells.csv.gz`` or ``*cells.parquet``).
    library_id : str or None
        Stored in ``uns["spatial"]`` for compatibility with scanpy / squidpy plotting.
    genes_only : bool
        Keep only features of type ``"Gene Expression"`` (drop control probes / codewords).
    verbose : bool

    Returns
    -------
    AnnData
        ``X`` = raw counts (sparse), ``obs`` = cell summary table (``cell_area``,
        ``transcript_counts`` ...), ``obs["imagecol"]`` / ``obs["imagerow"]`` and
        ``obsm["spatial"]`` = centroid in micrometres, ``var["gene_ids"]``, ``var["feature_types"]``.

    Examples
    --------
    >>> adata = sk.pp.read_xenium("Xenium_FFPE_Human_Breast_Cancer_Rep1/")
    >>> adata = sk.pp.read_xenium("Rep1_cell_feature_matrix.h5", cells_file="Rep1_cells.csv.gz")
    """
    if os.path.isdir(path):
        h5s = sorted(glob.glob(os.path.join(path, "*cell_feature_matrix.h5")))
        if not h5s:
            raise FileNotFoundError(f"no *cell_feature_matrix.h5 in {path}")
        h5 = h5s[0]
        if cells_file is None:
            cands = sorted(glob.glob(os.path.join(path, "*cells.parquet"))) + \
                    sorted(glob.glob(os.path.join(path, "*cells.csv.gz"))) + \
                    sorted(glob.glob(os.path.join(path, "*cells.csv")))
            if not cands:
                raise FileNotFoundError(f"no *cells.parquet / *cells.csv(.gz) in {path}")
            cells_file = cands[0]
    else:
        h5 = path
        if cells_file is None:
            raise ValueError("cells_file is required when path is the h5 file")

    X, barcodes, var = _read_10x_h5(h5)
    if cells_file.endswith(".parquet"):
        cells = pd.read_parquet(cells_file)
    else:
        cells = pd.read_csv(cells_file)
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    missing = set(barcodes) - set(cells.index)
    if missing:
        raise ValueError(f"{len(missing)} barcodes of the matrix are missing from {cells_file}")
    cells = cells.loc[barcodes]

    adata = ad.AnnData(X, obs=cells, var=var)
    adata.var_names_make_unique()
    if genes_only and "feature_types" in adata.var:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    adata.obs["imagecol"] = adata.obs["x_centroid"].astype(float).to_numpy()
    adata.obs["imagerow"] = adata.obs["y_centroid"].astype(float).to_numpy()
    adata.obsm["spatial"] = adata.obs[["imagecol", "imagerow"]].to_numpy()
    lib = library_id or os.path.basename(os.path.dirname(os.path.abspath(h5))) or "Xenium"
    adata.uns["spatial"] = {lib: {"scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 10.0},
                                  "metadata": {"source": "skny.pp.read_xenium", "coordinate_unit": "um"}}}
    if verbose:
        print(f"read_xenium: {adata.n_obs} cells x {adata.n_vars} genes from {os.path.basename(h5)}")
    return adata
