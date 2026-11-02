"""Stacked bar plot of cell density per sample for core / margin / dist."""
import numpy as np
import matplotlib.pyplot as plt

from ..preprocessing.calculate_region import REGION_ORDER


def region_density_plot(
    df, groups=None, colors=None, regions=None, samples=None, min_area_mm2=0.5,
    ylabel="density (cells / mm$^2$)", title=None, figsize=None, show=True,
):
    """Stacked bars of density per sample, one panel per region (core / margin / dist).

    Regions whose area is below ``min_area_mm2`` are not shown (their densities are
    unstable); the column is covered with a translucent hatch instead.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`skny.pp.region_density`.
    groups : list of str or None
        Groups to stack (order = stacking order). Default: all except ``_all_cells``.
    colors : dict or None
        Group -> colour. Default: matplotlib ``tab20``.
    regions : list of str or None
        Default: the regions present in ``df`` in core / margin / dist order.
    samples : list of str or None
        Sample order. Default: sorted.
    min_area_mm2 : float
        Hatch threshold (set to 0 to always show values).

    Returns
    -------
    (matplotlib Figure, array of Axes)
    """
    if groups is None:
        groups = [g for g in df["group"].unique() if g != "_all_cells"]
    if regions is None:
        regions = [r for r in REGION_ORDER if r in set(df["region"])]
    if samples is None:
        samples = sorted(df["sample"].unique())
    if colors is None:
        cmap = plt.get_cmap("tab20")
        colors = {g: cmap(i % 20) for i, g in enumerate(groups)}
    sub = df[df["group"].isin(groups)]
    d = sub.pivot_table(index=["sample", "region"], columns="group", values="density", fill_value=0.0, aggfunc="sum")
    d = d.reindex(columns=[g for g in groups if g in d.columns], fill_value=0.0)
    area = df.drop_duplicates(["sample", "region"]).set_index(["sample", "region"])["area_mm2"]
    small = {(smp, r): float(area.get((smp, r), 0.0)) < min_area_mm2 for smp in samples for r in regions}
    valid_tot = [float(d.loc[(smp, r)].sum()) for smp in samples for r in regions if (smp, r) in d.index and not small[(smp, r)]]
    ymax = (max(valid_tot) if valid_tot else 1.0) * 1.08 or 1.0
    if figsize is None:
        figsize = (max(6, 0.5 * len(samples) + 3), 3.5 * len(regions))
    fig, axes = plt.subplots(len(regions), 1, figsize=figsize, sharex=True, squeeze=False)
    axes = axes[:, 0]
    xs = np.arange(len(samples))
    for ax, r in zip(axes, regions):
        dd = d[d.index.get_level_values("region") == r].droplevel("region").reindex(samples).fillna(0.0)
        for smp in samples:
            if small[(smp, r)]:
                dd.loc[smp] = 0.0
        bottom = np.zeros(len(samples))
        for g in dd.columns:
            ax.bar(xs, dd[g].values, 0.8, bottom=bottom, color=colors.get(g), label=g if r == regions[0] else None)
            bottom += dd[g].values
        for k, smp in enumerate(samples):
            if small[(smp, r)]:
                ax.bar(k, ymax, 0.8, color="gray", alpha=0.3, hatch="//", edgecolor="gray", linewidth=0)
        ax.set_ylim(0, ymax); ax.set_ylabel(ylabel); ax.set_title(r, loc="left", fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    axes[0].legend(bbox_to_anchor=(1.01, 1), fontsize=8, loc="upper left")
    axes[-1].set_xticks(xs); axes[-1].set_xticklabels(samples, rotation=90, fontsize=8)
    if title is None:
        title = f"hatched column: region area < {min_area_mm2} mm$^2$ (not shown)"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    if show:
        plt.show()
    return fig, axes
