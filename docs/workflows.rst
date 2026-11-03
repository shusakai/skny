============================
Which workflow should I use?
============================

SKNY has two independent workflows. Both start from an ``AnnData`` of a tissue section and both
measure "how far from the tumor" something is, but they differ in the unit of analysis, in how the
distance is defined and in what comes out.

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * -
     - Workflow 1: tumor contour & distance bands
     - Workflow 2: tumor core / margin / dist
   * - Main function
     - :func:`skny.pp.calculate_distance`
     - :func:`skny.pp.calculate_region`
   * - Input
     - 10 µm grid built with :func:`skny.pp.grid` (counts summed per bin)
     - cell-level ``AnnData`` (one row per cell), coordinates in micrometres
   * - Tumor definition
     - marker genes (``pos_marker_ls`` / ``neg_marker_ls``) **or** cell labels
       (``use_label`` + ``pos_label_ls``)
     - cell label only (``ref_label``)
   * - Tumor region
     - binary tumor image smoothed with ``filters`` / ``threshold`` and contoured with OpenCV:
       a *spatial domain* with an explicit outline
     - bins within ``nest_radius`` (15 µm) of a tumor cell; no contour extraction, no hole filling
   * - Distance
     - shortest path from the contour on the grid; negative inside, positive outside;
       30 µm bands ``(-30, 0]``, ``(0, 30]``, ...
     - Euclidean distance to the nearest tumor cell; three regions ``core`` /
       ``margin`` (within ``margin`` = 100 µm of the core) / ``dist``
   * - Tissue mask
     - none (empty bins stay on the grid)
     - built from the density of all cells; empty space is ``outside``, never ``dist``
   * - Output
     - per bin: ``grid.shortest`` DataFrame, ``grid.uns[...]`` (images, contour, bands)
     - per cell: ``obs["region"]``, ``obs["region_dist_um"]``; per-sample maps and areas in
       ``uns["region"]``
   * - Downstream
     - :func:`skny.pl.distance_gene_plot`, :func:`skny.pl.clustering_heatmap`,
       :func:`skny.pl.deg_manhattan_plot`, :func:`skny.pp.convert_indivisual_solid` /
       :func:`skny.pp.convert_indivisual_peri` (single-TME objects per tumor nest)
     - :func:`skny.pp.region_density` (cells / mm²), :func:`skny.pl.region_map`,
       :func:`skny.pl.region_density_plot`
   * - Typical question
     - How does expression change with the distance from the tumor boundary? Which genes
       share a distance profile? How do individual tumor nests differ?
     - How many T cells / macrophages are in the core, the margin and the distant stroma,
       and how does that compare between samples?
   * - Run time (one Xenium section, ~170k cells)
     - about 1.5 min
     - seconds
   * - Reference
     - the method of the paper (Sakai *et al.*, PLOS Comput Biol 2025)
     - added in skny 0.3.0

Rules of thumb
--------------

* You have **marker genes** but no annotation, or you want **expression gradients, DEG or
  per-nest analysis**: workflow 1 (:doc:`workflow_contour`).
* You have an **annotated** object and want **cell densities per region**, across many
  samples: workflow 2 (:doc:`workflow_region`).
* You have an annotation and want both: pass the same label to both workflows. The tutorials
  :doc:`notebooks/region_analysis` (workflow 2) and :doc:`notebooks/label_based_contour`
  (workflow 1 from the same label) use the same Xenium section.

Minimal code
------------

.. code-block:: python

   import skny as sk

   # Workflow 1: grid -> contour -> distance bands
   adata = sk.pp.read_xenium("path/to/xenium_output")      # or any AnnData with cell coordinates
   grid = sk.pp.grid(adata, bin_um=10)
   grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"])
   sk.pl.distance_gene_plot(grid, "CD3E")

   # Workflow 2: cell label -> core / margin / dist -> density
   sk.pp.calculate_region(adata, "Tumor", label_key="cell_type", sample_key="sample")
   df = sk.pp.region_density(adata, "cell_type", sample_key="sample")
   sk.pl.region_density_plot(df, groups=["T cell", "Macrophage"])
