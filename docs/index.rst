SKNY (SpatialKNifeY) documentation
==================================

SKNY analyses the tumor microenvironment in single-cell spatial omics data (Xenium, CosMx,
PhenoCycler, any AnnData) by measuring **distance from the tumor**. It offers two workflows that
answer different questions; see :doc:`workflows` to choose.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Workflow 1: tumor contour & distance bands
     - Workflow 2: tumor core / margin / dist
   * - ``sk.pp.grid`` → ``sk.pp.calculate_distance``
     - ``sk.pp.calculate_region`` → ``sk.pp.region_density``
   * - 10 µm **grid**; tumor from **marker genes** or **cell labels**; smoothed image contour;
       shortest-path distance from the contour in 30 µm bands (negative inside, positive outside).
     - **Cell-level** data with coordinates in µm; tumor from a **cell label**; every cell gets
       ``core`` / ``margin`` / ``dist`` and its distance to the nearest tumor cell.
   * - Expression vs distance, gene clustering, DEG between bands, per-nest single-TME objects.
       The method of the paper.
     - Cell density (cells / mm²) per region and sample, region maps. Fast, no grid.
   * - :doc:`workflow_contour`
     - :doc:`workflow_region`

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   readme
   installation
   workflows

.. toctree::
   :maxdepth: 1
   :caption: Workflow 1 - Tumor contour (grid)

   workflow_contour
   notebooks/Xenium_analysis
   notebooks/label_based_contour
   notebooks/image_processing_options
   notebooks/single-TME_analysis
   notebooks/SKNY_paper_Sakai_et_al

.. toctree::
   :maxdepth: 1
   :caption: Workflow 2 - Tumor core / margin / dist (cells)

   workflow_region
   notebooks/region_analysis

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api
   history
   contributing
   authors

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
