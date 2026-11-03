"""Spatial expression map of one gene (cell-level or gridded AnnData), matplotlib."""
import numpy as np
import matplotlib.pyplot as plt


def gene_plot(
    adata, gene_symbols, ax=None, zoom_coord=None, vmin=None, vmax=None, cmap="viridis",
    show_color_bar=True, size=None, x_key="imagecol", y_key="imagerow", layer=None,
    title=None, figsize=(6, 6), show=False,
):
    """Draw the expression of a gene on the tissue coordinates (replacement for ``stlearn.pl.gene_plot``).

    Works for cell-level objects (e.g. from :func:`skny.pp.read_xenium`) and for gridded
    objects from :func:`skny.pp.grid` / ``stlearn.tl.cci.grid``. The y axis is inverted so that
    the orientation matches image coordinates.

    Parameters
    ----------
    adata : AnnData
        With ``obs[x_key]`` / ``obs[y_key]``.
    gene_symbols : str or list of str
        Gene(s); several genes are summed.
    ax : matplotlib Axes or None
    zoom_coord : list of 4 numbers or None
        ``[x_min, x_max, y_min, y_max]`` region of interest (same unit as the coordinates).
    vmin, vmax : float or None
        Colour limits.
    cmap : str
    show_color_bar : bool
    size : float or None
        Marker size; default is derived from the bin size for gridded data.
    layer : str or None
        Use ``adata.layers[layer]`` instead of ``X``.
    show : bool
        Call ``plt.show()``.

    Returns
    -------
    matplotlib Axes
    """
    genes = [gene_symbols] if isinstance(gene_symbols, str) else list(gene_symbols)
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"gene(s) not found: {missing}")
    M = adata[:, genes].layers[layer] if layer is not None else adata[:, genes].X
    vals = np.asarray(M.sum(axis=1)).ravel() if hasattr(M, "toarray") else np.asarray(M).sum(axis=1)
    x = adata.obs[x_key].to_numpy(dtype=float)
    y = adata.obs[y_key].to_numpy(dtype=float)

    if zoom_coord is not None:
        x0, x1, y0, y1 = zoom_coord
        m = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        x, y, vals = x[m], y[m], vals[m]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if size is None:
        if "grid_xedges" in adata.uns:
            bin_w = float(np.diff(adata.uns["grid_xedges"]).mean())
            # convert the bin width (data units) into points^2 for a square marker
            ax.set_xlim(x.min() - bin_w, x.max() + bin_w); ax.set_ylim(y.min() - bin_w, y.max() + bin_w)
            ax.figure.canvas.draw()
            px_per_unit = ax.transData.transform((bin_w, 0))[0] - ax.transData.transform((0, 0))[0]
            size = max((px_per_unit * 72.0 / ax.figure.dpi) ** 2, 1.0)
        else:
            size = 1.0
    order = np.argsort(vals)
    sc = ax.scatter(x[order], y[order], c=vals[order], s=size, marker="s", cmap=cmap, vmin=vmin, vmax=vmax,
                    linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title if title is not None else " + ".join(genes))
    if show_color_bar:
        ax.figure.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    if show:
        plt.show()
    return ax
