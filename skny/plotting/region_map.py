"""Plot the tumor core / margin / dist map of one sample."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from ..preprocessing.calculate_region import REGION_CODES, REGION_ORDER

REGION_COLORS = {"core": "#b2182b", "margin": "#fddbc7", "dist": "#4393c3", "outside": "#ffffff"}


def region_map(
    adata, sample=None, region_key="region", sample_key=None, spatial_key="spatial", coord_keys=None,
    highlight=None, highlight_key=None, highlight_color="lime", show_ref=True, colors=None,
    s=0.3, ax=None, figsize=(7, 7), title=None, show=True,
):
    """Draw the region map of one sample with the reference (tumor) cells on top.

    Parameters
    ----------
    adata : AnnData
        After :func:`skny.pp.calculate_region`.
    sample : str or None
        Sample to draw (required when regions were computed per sample).
    region_key : str
        ``key_added`` used in :func:`skny.pp.calculate_region`.
    sample_key, spatial_key, coord_keys
        As in :func:`skny.pp.calculate_region` (defaults are taken from ``uns[region_key]``).
    highlight : str, list of str or None
        Values of ``obs[highlight_key]`` (default: the label column used for the regions)
        to draw in ``highlight_color`` (e.g. immune cell types).
    show_ref : bool
        Draw the reference cells in black.
    colors : dict or None
        Region -> colour. Default :data:`REGION_COLORS`.
    ax : matplotlib Axes or None
    show : bool
        Call ``plt.show()``.

    Returns
    -------
    matplotlib Axes
    """
    info = adata.uns[region_key]
    params = info["params"]
    if sample_key is None:
        sample_key = params.get("sample_key")
    maps = info["maps"]
    if sample is None:
        if len(maps) != 1:
            raise ValueError(f"sample must be one of {list(maps)}")
        sample = next(iter(maps))
    if str(sample) not in maps:
        raise KeyError(f"no region map for sample {sample!r}; available: {list(maps)}")
    R = maps[str(sample)]
    colors = dict(REGION_COLORS, **(colors or {}))
    img = np.ones(tuple(np.asarray(R["region"]).shape) + (3,))
    for r, c in REGION_CODES.items():
        img[np.asarray(R["region"]) == c] = matplotlib.colors.to_rgb(colors[r])
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img, extent=np.asarray(R["extent"], float), origin="upper", interpolation="nearest")

    obs = adata.obs
    m = (obs[sample_key].astype(str) == str(sample)).to_numpy() if sample_key is not None else np.ones(adata.n_obs, bool)
    if coord_keys is None and spatial_key in adata.obsm:
        xy = np.asarray(adata.obsm[spatial_key], float)[:, :2]
    else:
        ck = coord_keys if coord_keys is not None else ("x_um", "y_um")
        xy = obs[list(ck)].to_numpy(float)
    label_key = params["label_key"]
    if show_ref:
        isref = obs[label_key].astype(str).isin(params["ref_label"]).to_numpy() & m
        ax.scatter(xy[isref, 0], xy[isref, 1], s=s / 2, c="k", alpha=0.35, linewidths=0, label="reference cells")
    if highlight is not None:
        hk = highlight_key or label_key
        hl = [highlight] if isinstance(highlight, str) else list(highlight)
        ish = obs[hk].astype(str).isin([str(v) for v in hl]).to_numpy() & m
        ax.scatter(xy[ish, 0], xy[ish, 1], s=s, c=highlight_color, alpha=0.6, linewidths=0, label=", ".join(map(str, hl)))
    area = R["area_mm2"]
    if title is None:
        title = f"{sample}: " + ", ".join(f"{r} {float(area[r]):.2f} mm$^2$" for r in REGION_ORDER if r in area)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(R["extent"][0], R["extent"][1]); ax.set_ylim(R["extent"][2], R["extent"][3])
    if show_ref or highlight is not None:
        ax.legend(loc="lower right", fontsize=7, markerscale=10)
    if show:
        plt.show()
    return ax
