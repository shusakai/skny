============================================
Workflow 1: tumor contour and distance bands
============================================

This is the method described in the SKNY paper. The tissue is divided into a 10 µm grid, the
tumor bins are turned into a binary image, smoothed and contoured, and the distance of every bin
from that contour is measured along the grid. Expression can then be studied as a function of the
distance from the tumor boundary, inside (negative distances) and outside (positive distances).

Steps
-----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - step
     - function
   * - 1. Read the data
     - :func:`skny.pp.read_xenium` (10x Xenium output), or any ``AnnData`` with per-cell
       coordinates in ``obs`` or ``obsm``
   * - 2. Grid
     - :func:`skny.pp.grid` sums the counts of the cells in each ``bin_um`` (10 µm) bin.
       ``use_label="cell_type"`` additionally stores the label fractions per bin.
   * - 3. Tumor contour and distance
     - :func:`skny.pp.calculate_distance`: tumor bins from ``pos_marker_ls`` (marker mode) or
       ``use_label`` + ``pos_label_ls`` (label mode); image smoothing with ``filters`` /
       ``threshold`` (:func:`skny.pp.apply_filters`); contour; shortest-path distance; 30 µm
       bands in ``grid.shortest``
   * - 4. Inspect
     - :func:`skny.pl.preprocessed_images`, :func:`skny.pl.equidistant_curves`,
       :func:`skny.pl.gene_plot`
   * - 5. Expression vs distance
     - :func:`skny.pl.distance_gene_plot` (one gene), :func:`skny.pl.clustering_heatmap`
       (all genes clustered by their distance profile), :func:`skny.pl.deg_manhattan_plot`
       (DEG between distance bands)
   * - 6. Single-TME analysis
     - :func:`skny.pp.convert_indivisual_solid` / :func:`skny.pp.convert_indivisual_peri`:
       one pseudo-bulk profile per tumor nest (inside the nest / band around it)

Two ways to define the tumor
----------------------------

.. code-block:: python

   grid = sk.pp.grid(adata, bin_um=10)
   grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"], neg_marker_ls=["KRT15"])   # marker mode

   grid = sk.pp.grid(adata, bin_um=10, use_label="cell_type")
   grid = sk.pp.calculate_distance(grid, use_label="cell_type", pos_label_ls=["Tumor"],
                                   neg_label_ls=["Myoepithelial"], min_fraction=0.5)          # label mode

Exactly one of the two must be given. The definition used is stored in
``grid.uns["tumor_definition"]``, the smoothing in ``grid.uns["image_processing"]``.

Tutorials
---------

* :doc:`notebooks/Xenium_analysis` - the full workflow on a 10x Xenium breast cancer section
  (marker mode)
* :doc:`notebooks/label_based_contour` - the same workflow starting from a cell-type annotation
  (label mode)
* :doc:`notebooks/image_processing_options` - choosing the smoothing of the tumor image
  (``filters`` / ``threshold``)
* :doc:`notebooks/single-TME_analysis` - per-nest analysis with ``convert_indivisual_*``
* :doc:`notebooks/SKNY_paper_Sakai_et_al` - reproduction of the paper figures

Grids built with stlearn (``st.tl.cci.grid``) are still accepted; bins are named ``grid_<n>``
with ``n = col * n_row + row``, numbered from 0 (use :func:`skny.pp.bin_names` to list them).
