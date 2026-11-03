=============
API reference
=============

Import SKNY as::

   import skny as sk

Preprocessing functions live in ``sk.pp``, plotting functions in ``sk.pl``. They are listed
below by workflow (see :doc:`workflows`).

Workflow 1: tumor contour and distance bands (grid)
===================================================

Reading and gridding
--------------------

.. autofunction:: skny.pp.read_xenium
.. autofunction:: skny.pp.grid
.. autofunction:: skny.pp.bin_names
.. autofunction:: skny.pp.grid_shape

Tumor contour and distance
--------------------------

.. autofunction:: skny.pp.calculate_distance
.. autofunction:: skny.pp.apply_filters

Single-TME objects (one per tumor nest)
---------------------------------------

.. autofunction:: skny.pp.convert_indivisual_solid
.. autofunction:: skny.pp.convert_indivisual_peri

Plotting (plotly)
-----------------

.. autofunction:: skny.pl.preprocessed_images
.. autofunction:: skny.pl.equidistant_curves
.. autofunction:: skny.pl.distance_gene_plot
.. autofunction:: skny.pl.clustering_heatmap
.. autofunction:: skny.pl.deg_manhattan_plot

Workflow 2: tumor core / margin / dist (cells)
==============================================

Regions and densities
---------------------

.. autofunction:: skny.pp.calculate_region
.. autofunction:: skny.pp.regions_from_coords
.. autofunction:: skny.pp.region_density

Plotting (matplotlib)
---------------------

.. autofunction:: skny.pl.region_map
.. autofunction:: skny.pl.region_density_plot

Shared
======

.. autofunction:: skny.pl.gene_plot
