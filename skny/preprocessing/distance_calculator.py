# load library
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import cv2
from PIL import Image
import networkx as nx
import math
from .grid import bin_names


# function
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
      

def cv2pil(image):
    ''' OpenCV type -> PIL type '''
    new_image = image.copy()
    if new_image.ndim == 2: 
        pass
    elif new_image.shape[2] == 3: 
        new_image = cv2.cvtColor(new_image, cv2.COLOR_BGR2RGB)
    elif new_image.shape[2] == 4: 
        new_image = cv2.cvtColor(new_image, cv2.COLOR_BGRA2RGBA)
    new_image = Image.fromarray(new_image)
    return new_image

DEFAULT_FILTERS = (("median", 3), ("box", 3))
_FILTER_OPS = ("median", "box", "mean", "blur", "gaussian", "open", "close", "erode", "dilate",
               "fill_holes", "min_area", "threshold")


def _parse_step(step):
    """Normalise a filter step to (op, args)."""
    if isinstance(step, str):
        return step.lower(), ()
    if isinstance(step, dict):
        step = dict(step); op = str(step.pop("op")).lower()
        return op, tuple(step.values())
    step = tuple(step)
    return str(step[0]).lower(), tuple(step[1:])


def apply_filters(mask, filters=DEFAULT_FILTERS, threshold=90):
    """Smooth a binary tumor mask with a sequence of image-processing steps and re-binarise it.

    The mask is handled as an 8-bit image (0 = background, 255 = tumor); linear filters blur it,
    morphological steps reshape it, and the final (or an explicit ``"threshold"``) step binarises
    it at ``threshold``. Used by :func:`skny.pp.calculate_distance`; call it directly to preview
    the effect of a filter set on ``grid.uns["marker"]``.

    Parameters
    ----------
    mask : ndarray (n_row, n_col)
        Boolean mask, or an 8-bit / colour image from ``grid.uns["marker"]`` (non-zero = tumor).
    filters : sequence of steps, or None
        Each step is ``"op"`` or ``("op", *args)``:

        * ``("median", k)``: median filter, ``k`` odd (default 3)
        * ``("box", k)`` / ``"mean"`` / ``"blur"``: box (mean) filter
        * ``("gaussian", k[, sigma])``: Gaussian filter (``sigma`` 0 = derived from ``k``)
        * ``("open", k[, iterations])``, ``("close", k[, iterations])``, ``("erode", k)``,
          ``("dilate", k)``: morphological operations with an elliptical kernel
        * ``"fill_holes"``: fill enclosed holes (binarises first)
        * ``("min_area", n)``: drop connected components smaller than ``n`` bins (binarises first)
        * ``"threshold"``: binarise at ``threshold`` at this point

        Default ``(("median", 3), ("box", 3))`` (the setting of the paper). ``None`` or ``()``
        returns the mask unchanged (useful to look at the raw mask; not accepted by
        :func:`skny.pp.calculate_distance`, where at least one smoothing step is required).
    threshold : int
        Value (0-255) above which a filtered pixel counts as tumor. With the default 3 x 3 box
        filter, 90 means at least 4 of the 9 pixels of the neighbourhood.

    Returns
    -------
    ndarray of bool (n_row, n_col)
        The filtered tumor mask.

    Examples
    --------
    >>> m = sk.pp.apply_filters(grid.uns["marker"], filters=[("gaussian", 5), ("close", 3), "fill_holes"], threshold=128)
    """
    m = np.asarray(mask)
    if m.ndim == 3:
        m = m.any(axis=2)
    img = np.where(m.astype(bool), 255, 0).astype(np.uint8)
    thr = int(threshold)

    def _binarise(a):
        return np.where(a > thr, 255, 0).astype(np.uint8)

    for step in (filters or ()):
        op, args = _parse_step(step)
        if op not in _FILTER_OPS:
            raise ValueError(f"unknown filter step {step!r}; allowed ops: {_FILTER_OPS}")
        k = int(args[0]) if args else 3
        if op == "threshold":
            img = _binarise(img)
        elif op == "median":
            if k % 2 == 0 or k < 3:
                raise ValueError("median kernel size must be odd and >= 3")
            img = cv2.medianBlur(img, k)
        elif op in ("box", "mean", "blur"):
            img = cv2.blur(img, (k, k))
        elif op == "gaussian":
            if k % 2 == 0 or k < 1:
                raise ValueError("gaussian kernel size must be odd and >= 1")
            sigma = float(args[1]) if len(args) > 1 else 0.0
            img = cv2.GaussianBlur(img, (k, k), sigma)
        elif op in ("open", "close", "erode", "dilate"):
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
            it = int(args[1]) if len(args) > 1 else 1
            cv_op = {"open": cv2.MORPH_OPEN, "close": cv2.MORPH_CLOSE, "erode": cv2.MORPH_ERODE, "dilate": cv2.MORPH_DILATE}[op]
            img = cv2.morphologyEx(img, cv_op, kernel, iterations=it)
        elif op == "fill_holes":
            from scipy import ndimage as ndi
            img = np.where(ndi.binary_fill_holes(img > thr), 255, 0).astype(np.uint8)
        elif op == "min_area":
            n_lab, lab, stats, _ = cv2.connectedComponentsWithStats(_binarise(img), connectivity=8)
            keep = np.zeros(n_lab, dtype=bool)
            keep[1:] = stats[1:, cv2.CC_STAT_AREA] >= k
            img = np.where(keep[lab], 255, 0).astype(np.uint8)
    return img > thr


def _tumor_bins(grid, pos_marker_ls=None, neg_marker_ls=None, use_label=None, pos_label_ls=None,
                neg_label_ls=None, min_cells=1, min_fraction=0.0):
    """Return a Series (index = grid.obs_names) with 1 for tumor bins and -1 otherwise."""
    if (pos_marker_ls is None) == (pos_label_ls is None):
        raise ValueError("give exactly one of pos_marker_ls (marker genes) or pos_label_ls (cell labels)")
    if isinstance(pos_marker_ls, str):
        pos_marker_ls = [pos_marker_ls]
    if isinstance(neg_marker_ls, str):
        neg_marker_ls = [neg_marker_ls]
    if isinstance(pos_label_ls, str):
        pos_label_ls = [pos_label_ls]
    if isinstance(neg_label_ls, str):
        neg_label_ls = [neg_label_ls]

    if pos_marker_ls is not None:
        # ---- expression-based: a bin is tumor if the summed expression of the positive markers is >= 1
        df = grid.to_df()
        missing = [g for g in list(pos_marker_ls) + list(neg_marker_ls or []) if g not in df.columns]
        if missing:
            raise KeyError(f"marker gene(s) not in grid.var_names: {missing}")
        pos = df[list(pos_marker_ls)].sum(axis=1).to_numpy(float) >= 1
        if neg_marker_ls:
            pos &= ~(df[list(neg_marker_ls)].sum(axis=1).to_numpy(float) >= 1)
        definition = dict(mode="marker", pos_marker_ls=list(pos_marker_ls), neg_marker_ls=list(neg_marker_ls or []))
    else:
        # ---- label-based: a bin is tumor if it contains >= min_cells cells of the positive labels
        #      (and at least min_fraction of its cells), and no (>= min_cells) cells of the negative labels
        if use_label is None:
            raise ValueError("use_label (the obs column given to skny.pp.grid(..., use_label=...)) is required with pos_label_ls")
        if use_label not in grid.uns:
            raise KeyError(
                f"uns['{use_label}'] not found: build the grid with skny.pp.grid(adata, use_label='{use_label}') "
                f"(or stlearn.tl.cci.grid(adata, ..., use_label='{use_label}'))"
            )
        if "n_cells" not in grid.obs:
            raise KeyError("obs['n_cells'] not found in the grid object")
        frac = pd.DataFrame(grid.uns[use_label]).loc[grid.obs_names]
        n_cells = grid.obs["n_cells"].to_numpy(float)
        missing = [l for l in list(pos_label_ls) + list(neg_label_ls or []) if l not in frac.columns]
        if missing:
            raise KeyError(f"label(s) not found in uns['{use_label}'] (available: {list(frac.columns)}): {missing}")
        f_pos = frac[list(pos_label_ls)].sum(axis=1).to_numpy(float)
        pos = (f_pos * n_cells >= min_cells) & (f_pos >= min_fraction)
        if neg_label_ls:
            f_neg = frac[list(neg_label_ls)].sum(axis=1).to_numpy(float)
            pos &= ~(f_neg * n_cells >= min_cells)
        definition = dict(mode="label", use_label=use_label, pos_label_ls=list(pos_label_ls),
                          neg_label_ls=list(neg_label_ls or []), min_cells=int(min_cells), min_fraction=float(min_fraction))
    series = pd.Series(np.where(pos, 1, -1), index=grid.obs_names, name="tumor_grid")
    return series, definition


def calculate_distance(
    grid, pos_marker_ls=None, neg_marker_ls=None,
    use_label=None, pos_label_ls=None, neg_label_ls=None, min_cells=1, min_fraction=0.0,
    filters=DEFAULT_FILTERS, threshold=90,
):
    """Detect the tumor spatial domain on a grid and compute the distance of every bin from its contour.

    This is the core step of the SKNY grid workflow (Sakai et al., PLOS Comput Biol 2025).
    The tumor bins are defined either from **marker gene expression** (``pos_marker_ls`` /
    ``neg_marker_ls``) or from **annotated cell labels** (``use_label`` + ``pos_label_ls`` /
    ``neg_label_ls``; the grid must have been built with ``skny.pp.grid(adata, use_label=...)``).
    The binary image is smoothed (``filters``: by default a 3 x 3 median filter followed by a
    3 x 3 box filter, binarised at ``threshold``), contoured with OpenCV, and the shortest-path
    distance from the contour to every bin is computed on a pixel graph (multi-source Dijkstra,
    8-connected). Distances are then binned into 30 um bands: negative bands lie inside the
    tumor, positive bands outside (peri-tumoral stroma).

    Parameters
    ----------
    grid : AnnData
        Gridded object from :func:`skny.pp.grid` or ``stlearn.tl.cci.grid`` (10 um bins
        recommended) with ``uns["grid_xedges"]`` and ``uns["grid_yedges"]``.
    pos_marker_ls : list of str or None
        Marker genes that define the tumor (e.g. ``["CDH1"]`` or ``["EPCAM", "KRT8"]``).
        A bin is tumor if the summed expression of these genes is >= 1.
    neg_marker_ls : list of str or None
        Genes whose expression (>= 1) excludes a bin from the tumor (e.g. stromal markers).
    use_label : str or None
        Name of the cell-label column given to ``skny.pp.grid(adata, use_label=...)``; the
        per-bin label fractions are read from ``grid.uns[use_label]``. Required with ``pos_label_ls``.
    pos_label_ls : list of str or None
        Cell labels that define the tumor (e.g. ``["Malignant"]`` or ``["DCIS", "IDC"]``).
        A bin is tumor if it contains at least ``min_cells`` such cells and they make up at least
        ``min_fraction`` of the cells in the bin.
    neg_label_ls : list of str or None
        Cell labels whose presence (>= ``min_cells`` cells) excludes a bin from the tumor.
    min_cells : int
        Minimum number of positive-label cells per bin (label mode only). Default 1.
    min_fraction : float
        Minimum fraction (0-1) of positive-label cells per bin (label mode only). Default 0,
        i.e. any bin touched by a tumor cell; use e.g. 0.5 for a stricter domain.
    filters : sequence of steps or None
        Image-processing steps applied to the binary tumor image before contouring, e.g.
        ``[("median", 3), ("box", 3)]`` (default), ``[("gaussian", 5), ("close", 3), "fill_holes"]``
        or ``[("median", 5), ("min_area", 10)]``. See :func:`skny.pp.apply_filters` for all
        steps. At least one smoothing step is required: without smoothing the contour breaks
        into thousands of fragments and the distance calculation becomes extremely slow.
    threshold : int
        Binarisation threshold (0-255) applied after the filters. Default 90.

    Returns
    -------
    AnnData
        The same ``grid`` with the processing images added to ``grid.uns``
        (``"marker"``, ``"marker_median"``, ``"marker_delineation"``,
        ``"shotest_30_delineation"``, ``"shotest_region_<band>_delineation"``), the tumor
        definition in ``grid.uns["tumor_definition"]``, the filter settings in
        ``grid.uns["image_processing"]``, and the per-bin distance table stored
        as the attribute ``grid.shortest`` (a DataFrame with a ``"region"`` column holding the
        distance band of each bin).

    See Also
    --------
    skny.pp.grid : builds the grid (pass ``use_label`` to enable the label mode).
    skny.pl.preprocessed_images, skny.pl.equidistant_curves, skny.pl.distance_gene_plot,
    skny.pl.clustering_heatmap, skny.pl.deg_manhattan_plot,
    skny.pp.convert_indivisual_solid, skny.pp.convert_indivisual_peri,
    skny.pp.calculate_region : cell-level core / margin / dist regions from a label (no grid).

    Examples
    --------
    Marker-based (as in the paper):

    >>> grid = sk.pp.grid(adata, bin_um=10)
    >>> grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"])

    Label-based, from a ``cell_type`` annotation in ``adata.obs``:

    >>> grid = sk.pp.grid(adata, bin_um=10, use_label="cell_type")
    >>> grid = sk.pp.calculate_distance(grid, use_label="cell_type", pos_label_ls=["Malignant"],
    ...                                 neg_label_ls=["Fibroblast"], min_fraction=0.5)
    >>> sk.pl.distance_gene_plot(grid, "KRT15")

    Different smoothing (Gaussian filter, closing, hole filling, stricter threshold):

    >>> grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"],
    ...                                 filters=[("gaussian", 5), ("close", 3), "fill_holes"], threshold=128)
    """
    ## get col and row length data from grid object --------------------------
    N_ROW = len(grid.uns["grid_yedges"]) - 1
    N_COL = len(grid.uns["grid_xedges"]) - 1
    
    ## Define tumor bins (marker genes or cell labels) ---------------------------
    pos_series, definition = _tumor_bins(
        grid, pos_marker_ls, neg_marker_ls, use_label, pos_label_ls, neg_label_ls, min_cells, min_fraction
    )
    grid.uns["tumor_definition"] = definition

    # expand information of all grid
    df_grid_tumor = pd.merge(
        pd.DataFrame(index=bin_names(grid)), 
        pos_series, right_index=True, left_index=True, how="left"
    ).fillna(-1)


    ## Generate the binary tumor image (one pixel per bin) ---------------------------
    # Tumor bins are drawn in the colour the original matplotlib/viridis rendering used
    # (RGB 253, 231, 36; stored as BGR for OpenCV), background in black. Building the image
    # directly avoids the interpolation / HiDPI artefacts of rendering through a figure.
    tumor_mask = np.asarray(df_grid_tumor["tumor_grid"], dtype=float).reshape(N_COL, N_ROW).T == 1
    img = np.zeros((N_ROW, N_COL, 3), dtype=np.uint8)
    img[tumor_mask] = (36, 231, 253)
    grid.uns["marker"] = img.copy() # add to grid object

    
    ## Smoothing (configurable filters) and binarisation ---------------------------
    filters = tuple(filters) if filters else ()
    if not filters:
        raise ValueError(
            "filters must contain at least one smoothing step, e.g. the default "
            "[('median', 3), ('box', 3)]; see skny.pp.apply_filters for the available steps"
        )
    smoothed = apply_filters(tumor_mask, filters=filters, threshold=threshold)
    # recorded as strings ("median 3") so that the object can be written to h5ad
    grid.uns["image_processing"] = dict(
        filters=[" ".join(str(x) for x in (op,) + tuple(args)) for op, args in map(_parse_step, filters)],
        threshold=int(threshold),
    )
    img_med = np.zeros((N_ROW, N_COL, 3), dtype=np.uint8)
    img_med[smoothed] = (36, 231, 253)
    grid.uns["marker_median"] = img_med.copy()  # smoothed tumor image (key kept for compatibility)

    ## Delineation --------------------------------
    img_binary = np.where(smoothed, 255, 0).astype(np.uint8)
    contours, hierarchy = cv2.findContours(img_binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contours_ = list(contours) if hierarchy is not None else []

    img_color_with_contours = cv2.drawContours(img_med, contours_, -1, (0,255,0), 1)
    grid.uns["marker_delineation"] = img_color_with_contours # add to grid object


    ## Converting image data to graph stracture ---------------------------------
    # OpenCV -> PIL
    image = cv2pil(img_color_with_contours)

    # extract size
    width, height = image.size
    # generate graph
    graph = nx.Graph()
    
    tumor_contour_pixel_ls = []
    tumor_pixel_ls = []
    edge_pixel_ls = []
    # ピクセルをノードに変換
    for y in range(height):
        for x in range(width):
            pixel_value = image.getpixel((x, y))
            graph.add_node((x, y), color=pixel_value)
    
            # extract tumor coutour coordination
            if pixel_value == (0, 255, 0):
                tumor_contour_pixel_ls += [(x, y)]
    
            # extract tumor coordination
            elif pixel_value == (253, 231, 36):
            #elif pixel_value == (255, 255, 0):

                tumor_pixel_ls += [(x, y)]
    
            if (y == 0) | (y == height-1) | (x == 0) | (x == width-1):
                edge_pixel_ls += [(x, y)]
    
    # get nodes
    nodes_ls = graph.nodes
    nodes_ls = list(nodes_ls)
    
    # check node inside contour
    inside_contour_ls = []
    for node in nodes_ls:
        for i in contours_:
            if cv2.pointPolygonTest(i, node, False) == 1:
                inside_contour_ls += [node]
                break
    
    # extract tumor inside contour
    tumor_pixel_ls = [i for i in tumor_pixel_ls if i in inside_contour_ls]
    
    # extract non-tumor inside contour
    #nontumor_pixel_inside_contour_ls = [i for i in inside_contour_ls if i not in tumor_pixel_ls]
    
    # 隣接ピクセル間にエッジを作成
    for y in range(height):
        for x in range(width):
            current_node = (x, y)
            neighbors = [
                (x, y - 1),  # 上
                (x, y + 1),  # 下
                (x - 1, y),  # 左
                (x + 1, y)   # 右
            ]
            # 上下左右のエッジ
            for neighbor in neighbors:
                if neighbor in graph.nodes:
                    graph.add_edge(current_node, neighbor)
    
            # 斜めのエッジ
            if (y != height) & (x != width):
                node1 = (x, y)
                node2 = (x + 1, y + 1)
    
                # ノードの座標からユークリッド距離を計算
                distance = math.sqrt((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)
                graph.add_edge(node1, node2, weight=distance)
    
                node1 = (x + 1, y)
                node2 = (x, y + 1)
    
                # ノードの座標からユークリッド距離を計算
                distance = math.sqrt((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)
                graph.add_edge(node1, node2, weight=distance)
    
    
    #remove non-tumor inside contour from graph
    #graph.remove_nodes_from(nontumor_pixel_inside_contour_ls)
    
    # delete contour at the edge of image 
    tumor_edge_contour_pixel_ls = list(set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    
    tumor_pixel_ls = list(
        set(tumor_pixel_ls) | (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    )
    
    tumor_contour_pixel_ls = list(
        set(tumor_contour_pixel_ls) - (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    )

    #tumor_contour_pixel_ls = [i for i in tumor_contour_pixel_ls if i in list(graph.nodes)]
    
    ## Calculate distance of each nodes from tumor surface ---------------------------------------------------
    shortest_paths = nx.multi_source_dijkstra(graph, tumor_contour_pixel_ls, cutoff=16, weight="weight")
    
    # shotest length from contour
    df_shotest = pd.DataFrame.from_dict(shortest_paths[0], orient="index", columns=["euclidean"])
    
    # negative value for inside contour
    df_shotest.loc[df_shotest.index[df_shotest.index.isin(tumor_pixel_ls)]] = df_shotest.loc[df_shotest.index[df_shotest.index.isin(tumor_pixel_ls)]] * -1
    df_shotest["euclidean_round"] = df_shotest["euclidean"].round(0)
    df_nodes = pd.DataFrame(index=nodes_ls)
    
    df_shotest = pd.merge(
        df_nodes, df_shotest, right_index=True, left_index=True, how="left"
    ).fillna(np.nan)

    ## Delineation ----------------------------------
    col_ls = []
    for i in df_shotest["euclidean_round"]:
        if i == 0:
            col_ls += [(0, 255, 0)]
        elif i < 0:
            col_ls += [(0, 255, 255)]
        else:
            col_ls += [(0, 0, 0)]
    
    col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
    grid.uns["marker_median_delineation"] = col_arr.copy() # add to grid object

    ## Delineation with 30μm interval -------------
    col_ls = []
    for i in df_shotest["euclidean_round"]:
        if i == 0:
            col_ls += [(0, 255, 0)]
        elif (i % 3 == 0) & (i > 0):
            col_ls += [(0, 0, 255)]
        elif (i % 3 == 0) & (i < 0):
            col_ls += [(0, 150, 255)]
        elif i < 0:
            col_ls += [(0, 255, 255)]
        else:
            col_ls += [(0, 0, 0)]
    
    col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
    grid.uns["shotest_30_delineation"] = col_arr.copy() # add to grid object

    ## Discretization ---------------------------------
    df_shotest["region"] = pd.cut(
        df_shotest.dropna()["euclidean"], 
        bins=list(range(-15, 16, 3))
    )
    
    for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index:
            
        col_ls = []
        for s in df_shotest["region"]:
            if s == i:
                col_ls += [(255, 255, 255)]
            else:
                col_ls += [(0, 0, 0)]
        
        col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
        grid.uns[f"shotest_region_{i}_delineation"] = col_arr # add to grid object
    
    # annotation of dataframe of shortest distance from surface
    #grid.uns["shortest"] = df_shotest
    setattr(grid, 'shortest', df_shotest)

    return grid
