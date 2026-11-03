======================================
Workflow 2: tumor core / margin / dist
======================================

Since v0.3.0 SKNY provides a second, contour-free way to define tumor regions.
:func:`skny.pp.calculate_region` works on **cell-level** data (one row per cell, e.g. a scanpy
object) whose cells have already been annotated, and defines the tumor from a **cell label**
such as ``"Malignant"``. Every cell is assigned to ``core`` / ``margin`` / ``dist`` and
:func:`skny.pp.region_density` turns that into cell densities (cells / mm²) per region and
sample. Compare with workflow 1 (:doc:`workflow_contour`), which contours the tumor on a
10 µm grid and reports expression as a function of the distance from the contour; see
:doc:`workflows` for when to use which.

Steps
-----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - step
     - function
   * - 1. Annotated cells
     - any ``AnnData`` with a cell-type column in ``obs`` and coordinates in micrometres
       (``obsm["spatial"]`` or two ``obs`` columns)
   * - 2. Regions
     - :func:`skny.pp.calculate_region` (per sample with ``sample_key``);
       :func:`skny.pp.regions_from_coords` for plain coordinate arrays
   * - 3. Densities
     - :func:`skny.pp.region_density`: cells / mm² of each group in each region
   * - 4. Plots
     - :func:`skny.pl.region_map`, :func:`skny.pl.region_density_plot`

Definitions
-----------

All computations use a ``bin_um`` (10 µm) grid per sample.

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - region
     - definition
   * - ``core``
     - bins within ``nest_radius`` (15 µm) of a reference (tumor) cell, i.e. the area actually
       occupied by tumor cells. No morphological closing or hole filling is applied, so stroma
       enclosed by tumor is **not** core. Connected components containing fewer than
       ``min_nest_cells`` (5) reference cells are discarded as isolated, probably mislabelled cells.
   * - ``margin``
     - tissue outside the core and within ``margin`` (100 µm) of the core boundary
       (peri-tumoral stroma)
   * - ``dist``
     - tissue farther than ``margin`` from the core boundary (distant stroma)
   * - ``outside``
     - outside the tissue mask; never counted

The **tissue mask** is the density of *all* cells smoothed with a Gaussian of
``tissue_smooth`` (60 µm) and thresholded at ``tissue_min_density`` (250 cells / mm²), holes
filled, fragments smaller than ``min_tissue_area`` (0.1 mm²) removed. Empty space (outside the
section, tears, the edge of a TMA core) is therefore never counted as ``dist``.

Distances are Euclidean distances to the nearest reference cell (``scipy.ndimage.distance_transform_edt``).
With ``margin_from="reference"`` the margin is measured from *any* reference cell instead of the
retained core boundary (legacy behaviour: isolated reference cells that were discarded from the
core still create a margin halo).

Minimal example
---------------

.. code-block:: python

   import skny as sk

   # adata: cells x genes, coordinates in micrometres in adata.obsm["spatial"],
   # cell types in adata.obs["cell_type"], sample id in adata.obs["sample"]
   sk.pp.calculate_region(adata, "Malignant", label_key="cell_type", sample_key="sample")

   adata.obs["region"]            # core / margin / dist / outside / not_computed
   adata.obs["region_dist_um"]    # distance to the nearest tumor cell
   adata.uns["region"]["maps"]    # per-sample region map, extent and areas (mm2)

   df = sk.pp.region_density(adata, "immune_subtype", sample_key="sample")   # cells / mm2
   sk.pl.region_map(adata, sample="S1", highlight=["T cell", "Macrophage"])
   sk.pl.region_density_plot(df, groups=["CD8+ T", "Treg", "M2-like"], min_area_mm2=0.5)

Coordinates may also be given as two ``obs`` columns: ``coord_keys=("x_um", "y_um")``.
Use ``key_added`` to keep several definitions side by side (e.g. ``key_added="region_morula"``
with ``core_only=True`` for a rare tumor subpopulation, compared with the main tumor core of the
same sample).

Output of ``region_density``
----------------------------

One row per sample × region × group with ``n`` (cells), ``area_mm2`` and ``density`` (= n / area).
The group ``_all_cells`` gives the total cell density of the region. Regions whose area is below
``min_area_mm2`` are hatched in :func:`skny.pl.region_density_plot` instead of showing an unstable
density.

Tutorial
--------

* :doc:`notebooks/region_analysis` - complete example on the 10x Xenium breast cancer dataset
  (annotation with scanpy, regions, densities per region)

The same annotation can also drive workflow 1 (:doc:`notebooks/label_based_contour`).
See :doc:`api` for all parameters.
