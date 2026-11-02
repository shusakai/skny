=============
API reference
=============

Import SKNY as::

   import skny as sk

Preprocessing: ``sk.pp``
------------------------

Expression-based tumor contour (grid data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skny.pp.calculate_distance
.. autofunction:: skny.pp.convert_indivisual_solid
.. autofunction:: skny.pp.convert_indivisual_peri

Cell-label based regions: tumor core / margin / dist (cell-level data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skny.pp.calculate_region
.. autofunction:: skny.pp.region_density
.. autofunction:: skny.pp.regions_from_coords

Plotting: ``sk.pl``
-------------------

.. autofunction:: skny.pl.preprocessed_images
.. autofunction:: skny.pl.equidistant_curves
.. autofunction:: skny.pl.clustering_heatmap
.. autofunction:: skny.pl.distance_gene_plot
.. autofunction:: skny.pl.deg_manhattan_plot
.. autofunction:: skny.pl.region_map
.. autofunction:: skny.pl.region_density_plot
