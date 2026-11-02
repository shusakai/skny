"""Cell density (cells / mm^2) of cell groups per tumor region."""
import numpy as np
import pandas as pd

from .calculate_region import REGION_ORDER


def region_density(adata, group_key, region_key="region", sample_key=None, groups=None, regions=None, include_all=True):
    """Count cells of each group per region and divide by the region area.

    Requires :func:`skny.pp.calculate_region` to have been run with ``key_added=region_key``.

    Parameters
    ----------
    adata : AnnData
        Cell-level object with ``obs[region_key]`` and ``uns[region_key]``.
    group_key : str
        ``obs`` column with the cell groups to count (e.g. immune subtypes).
        Cells with a missing value in this column are only counted in ``_all_cells``.
    region_key : str
        ``key_added`` used in :func:`skny.pp.calculate_region`.
    sample_key : str or None
        Sample column. Must match the one used in :func:`skny.pp.calculate_region`.
    groups : list of str or None
        Restrict / order the groups. Default: every value present.
    regions : list of str or None
        Default: ``["core", "margin", "dist"]`` (only ``core`` when ``core_only`` was used).
    include_all : bool
        Also report ``_all_cells`` (all cells in the region) per sample x region.

    Returns
    -------
    pandas.DataFrame
        Columns ``sample, region, group, n, area_mm2, density`` (``density`` = n / area, cells per mm^2;
        NaN when the area is 0). Every sample x region x group combination is present (n = 0 when absent).
    """
    if region_key not in adata.uns or region_key not in adata.obs:
        raise KeyError(f"run skny.pp.calculate_region(..., key_added='{region_key}') first")
    info = adata.uns[region_key]
    if sample_key is None:
        sample_key = info["params"].get("sample_key")
    if regions is None:
        regions = ["core"] if info["params"].get("core_only") else list(REGION_ORDER)
    obs = adata.obs
    samples = obs[sample_key].astype(str) if sample_key is not None else pd.Series("sample", index=obs.index)
    grp = obs[group_key]
    if groups is None:
        groups = [g for g in (grp.cat.categories if hasattr(grp, "cat") else pd.unique(grp.dropna()))]
    groups = [str(g) for g in groups]
    grp = grp.astype(object).where(grp.notna(), None).astype(str)
    reg = obs[region_key].astype(str)
    rows = []
    for s, R in info["maps"].items():
        area = R["area_mm2"]
        ms = (samples == s).to_numpy()
        for r in regions:
            a = float(area.get(r, 0.0)) if hasattr(area, "get") else float(area[r])
            mr = ms & (reg == r).to_numpy()
            cnt = grp[mr].value_counts()
            for g in groups:
                n = int(cnt.get(g, 0))
                rows.append(dict(sample=s, region=r, group=g, n=n, area_mm2=a, density=n / a if a > 0 else np.nan))
            if include_all:
                n = int(mr.sum())
                rows.append(dict(sample=s, region=r, group="_all_cells", n=n, area_mm2=a, density=n / a if a > 0 else np.nan))
    return pd.DataFrame(rows, columns=["sample", "region", "group", "n", "area_mm2", "density"])
