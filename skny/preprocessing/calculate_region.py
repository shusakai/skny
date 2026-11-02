"""Cell-label based tumor region analysis (core / margin / dist).

This is a variant of :func:`skny.pp.calculate_distance`.  Instead of
contouring the tumor from marker *gene expression* on a grid, the tumor is
defined from a *cell label* (e.g. ``cell_type == "Malignant"``) and the
regions are built from the Euclidean distance to the nearest reference cell:

* ``core``   : grid bins within ``nest_radius`` um of a reference cell
               (the area actually occupied by tumor cells; no closing or
               hole filling, so stroma enclosed by tumor is *not* core)
* ``margin`` : tissue outside the core but within ``margin`` um of the core
               boundary (peri-tumoral stroma)
* ``dist``   : tissue farther than ``margin`` um from the core boundary
* ``outside``: bins outside the tissue mask (not analysed)

The tissue mask is derived from the density of *all* cells so that empty
areas (outside the section, holes, TMA core edges) are never counted as
``dist``.  Regions are computed per sample when ``sample_key`` is given.
"""
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

REGION_ORDER = ["core", "margin", "dist"]
REGION_CODES = {"outside": 0, "core": 1, "margin": 2, "dist": 3}
NOT_COMPUTED = "not_computed"


def _remove_small(mask, min_px):
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = ndi.sum(mask, lab, range(1, n + 1))
    keep = np.zeros(n + 1, bool)
    keep[1:] = np.asarray(sizes) >= min_px
    return keep[lab]


def regions_from_coords(
    x, y, is_ref,
    bin_um=10.0, nest_radius=15.0, margin=100.0, min_nest_cells=5, min_ref_cells=20,
    tissue_smooth=60.0, tissue_min_density=250.0, min_tissue_area=0.1, core_only=False,
    margin_from="nest",
):
    """Compute the region map for one sample from cell coordinates.

    Parameters
    ----------
    x, y : array-like, shape (n_cells,)
        Cell coordinates in micrometres.
    is_ref : array-like of bool, shape (n_cells,)
        True for reference (tumor) cells.
    bin_um : float
        Grid size in um.
    nest_radius : float
        Bins within this distance (um) of a reference cell are ``core``.
    margin : float
        Width (um) of the ``margin`` band outside the core.
    min_nest_cells : int
        Connected core components containing fewer reference cells than this
        are discarded (isolated, probably mislabelled cells in the stroma).
    min_ref_cells : int
        If the sample has fewer reference cells than this, ``None`` is returned.
    tissue_smooth : float
        Sigma (um) of the Gaussian filter applied to the all-cell density
        image to build the tissue mask.
    tissue_min_density : float
        Minimum smoothed density (cells / mm^2) to be called tissue.
    min_tissue_area : float
        Tissue fragments smaller than this (mm^2) are removed from the mask.
    core_only : bool
        If True, only ``core`` is defined (margin / dist stay ``outside``).
    margin_from : {"nest", "reference"}
        ``"nest"`` (default): the margin band is measured from the boundary of the
        retained core (``0 < d_nest <= margin``). ``"reference"``: legacy behaviour,
        measured from *any* reference cell (``d_ref <= nest_radius + margin``), so
        isolated reference cells that were discarded from the core still create a
        margin halo around them.

    Returns
    -------
    dict or None
        ``region`` (2-D int8 array of :data:`REGION_CODES`), ``dist_um``
        (2-D float array: distance to the nearest reference cell),
        ``cell_region`` (array of str per cell), ``cell_dist_um`` (per cell),
        ``area_mm2`` (dict region -> mm^2), ``extent`` (x0, x1, y1, y0; for
        ``imshow(origin="upper")``), ``bin_um``.
    """
    x = np.asarray(x, float); y = np.asarray(y, float); is_ref = np.asarray(is_ref, bool)
    if len(x) == 0 or is_ref.sum() < min_ref_cells:
        return None
    pad = int(np.ceil((margin + nest_radius) / bin_um)) + 3
    x0, y0 = x.min() - pad * bin_um, y.min() - pad * bin_um
    nx = int((x.max() - x0) / bin_um) + pad + 1
    ny = int((y.max() - y0) / bin_um) + pad + 1
    ix = ((x - x0) / bin_um).astype(int); iy = ((y - y0) / bin_um).astype(int)

    # tissue mask from all-cell density
    allc = np.zeros((ny, nx)); np.add.at(allc, (iy, ix), 1)
    thr = tissue_min_density * (bin_um * bin_um) / 1e6          # cells per bin
    tissue = ndi.gaussian_filter(allc, tissue_smooth / bin_um) > thr
    tissue = ndi.binary_fill_holes(tissue)
    tissue = _remove_small(tissue, min_tissue_area * 1e6 / (bin_um * bin_um))

    # distance to nearest reference cell
    occ = np.zeros((ny, nx), bool); occ[iy[is_ref], ix[is_ref]] = True
    d_ref = ndi.distance_transform_edt(~occ) * bin_um
    nest = (d_ref <= nest_radius) & tissue
    lab, n = ndi.label(nest)
    if n:
        cnt = np.bincount(lab[iy[is_ref], ix[is_ref]], minlength=n + 1)
        keep = cnt >= min_nest_cells; keep[0] = False
        nest = keep[lab]
    if nest.sum() == 0:
        return None

    reg = np.zeros((ny, nx), np.int8)
    reg[nest] = REGION_CODES["core"]
    if not core_only:
        if margin_from == "nest":
            d_nest = ndi.distance_transform_edt(~nest) * bin_um    # distance to the core boundary
            near = d_nest <= margin
        elif margin_from == "reference":
            near = d_ref <= nest_radius + margin
        else:
            raise ValueError("margin_from must be 'nest' or 'reference'")
        reg[tissue & ~nest & near] = REGION_CODES["margin"]
        reg[tissue & ~near] = REGION_CODES["dist"]
    names = np.array(list(REGION_CODES.keys()))
    area = {r: float((reg == c).sum()) * bin_um * bin_um / 1e6 for r, c in REGION_CODES.items() if r != "outside"}
    return dict(
        region=reg, dist_um=d_ref, cell_region=names[reg[iy, ix]], cell_dist_um=d_ref[iy, ix],
        area_mm2=area, extent=np.array([x0, x0 + nx * bin_um, y0 + ny * bin_um, y0]), bin_um=float(bin_um),
    )


def calculate_region(
    adata, ref_label, label_key="cell_type", sample_key=None, spatial_key="spatial", coord_keys=None,
    bin_um=10.0, nest_radius=15.0, margin=100.0, min_nest_cells=5, min_ref_cells=20,
    tissue_smooth=60.0, tissue_min_density=250.0, min_tissue_area=0.1, core_only=False,
    margin_from="nest", key_added="region", copy=False, verbose=True,
):
    """Assign every cell to tumor ``core`` / ``margin`` / ``dist`` from a cell label.

    Parameters
    ----------
    adata : AnnData
        Cell-level object (one row per cell) with coordinates in **micrometres**.
    ref_label : str or list of str
        Value(s) of ``adata.obs[label_key]`` that define the tumor (e.g. ``"Malignant"``).
    label_key : str
        ``obs`` column holding the cell labels.
    sample_key : str or None
        ``obs`` column identifying samples / sections. Regions are computed per sample.
    spatial_key : str
        Key of ``adata.obsm`` holding (x, y) coordinates in um (used when ``coord_keys`` is None).
    coord_keys : tuple of str or None
        Alternatively, two ``obs`` columns with x and y in um.
    bin_um, nest_radius, margin, min_nest_cells, min_ref_cells, tissue_smooth,
    tissue_min_density, min_tissue_area, core_only, margin_from
        See :func:`regions_from_coords`.
    key_added : str
        ``obs[key_added]`` receives the region label (``core`` / ``margin`` / ``dist`` /
        ``outside`` / ``not_computed``), ``obs[key_added + "_dist_um"]`` the distance to the
        nearest reference cell, ``uns[key_added]`` the parameters and per-sample maps.
    copy : bool
        Return a copy instead of modifying ``adata`` in place.

    Returns
    -------
    AnnData or None
        ``adata`` (if ``copy``), otherwise None. Results are stored in ``obs`` and ``uns``.
    """
    adata = adata.copy() if copy else adata
    if isinstance(ref_label, str):
        ref_label = [ref_label]
    if label_key not in adata.obs:
        raise KeyError(f"obs['{label_key}'] not found")
    if coord_keys is not None:
        xy = adata.obs[list(coord_keys)].to_numpy(float)
    else:
        if spatial_key not in adata.obsm:
            raise KeyError(f"obsm['{spatial_key}'] not found; pass coord_keys=(x_col, y_col) instead")
        xy = np.asarray(adata.obsm[spatial_key], float)[:, :2]
    is_ref = adata.obs[label_key].astype(str).isin([str(v) for v in ref_label]).to_numpy()
    samples = adata.obs[sample_key].astype(str).to_numpy() if sample_key is not None else np.array(["sample"] * adata.n_obs)
    params = dict(
        ref_label=list(map(str, ref_label)), label_key=label_key, sample_key=sample_key, bin_um=float(bin_um),
        nest_radius=float(nest_radius), margin=float(margin), min_nest_cells=int(min_nest_cells),
        min_ref_cells=int(min_ref_cells), tissue_smooth=float(tissue_smooth),
        tissue_min_density=float(tissue_min_density), min_tissue_area=float(min_tissue_area), core_only=bool(core_only),
        margin_from=str(margin_from),
    )
    cell_region = np.full(adata.n_obs, NOT_COMPUTED, dtype=object)
    cell_dist = np.full(adata.n_obs, np.nan)
    maps = {}
    for s in pd.unique(samples):
        m = samples == s
        R = regions_from_coords(
            xy[m, 0], xy[m, 1], is_ref[m], bin_um=bin_um, nest_radius=nest_radius, margin=margin,
            min_nest_cells=min_nest_cells, min_ref_cells=min_ref_cells, tissue_smooth=tissue_smooth,
            tissue_min_density=tissue_min_density, min_tissue_area=min_tissue_area, core_only=core_only,
            margin_from=margin_from,
        )
        if R is None:
            if verbose:
                print(f"[skny] {s}: skipped (fewer than {min_ref_cells} reference cells or no tumor nest)")
            continue
        cell_region[m] = R["cell_region"]; cell_dist[m] = R["cell_dist_um"]
        maps[str(s)] = dict(region=R["region"], extent=R["extent"], area_mm2=R["area_mm2"], bin_um=R["bin_um"])
        if verbose:
            print(f"[skny] {s}: " + ", ".join(f"{r}={R['area_mm2'][r]:.3f} mm2" for r in REGION_ORDER))
    cats = REGION_ORDER + ["outside", NOT_COMPUTED]
    adata.obs[key_added] = pd.Categorical(cell_region, categories=cats)
    adata.obs[key_added + "_dist_um"] = cell_dist
    adata.uns[key_added] = dict(params=params, codes=dict(REGION_CODES), maps=maps)
    return adata if copy else None
