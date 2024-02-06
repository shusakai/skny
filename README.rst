.. image:: https://img.shields.io/pypi/v/skny.svg
        :target: https://pypi.python.org/pypi/skny

.. image:: https://readthedocs.org/projects/skny/badge/?version=latest
        :target: https://skny.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

SKNY - SpatialKNifeY
=====================

**SKNY** is a tools for spatial anlysis stratified by distance from tumor for multiple platform such as Xenium, CosMx, and PhenoCycler. 
It automatically contours the tumor based on the spatial omics expression data and calculates the distance from this contour to each coordinate in space.
Using this distance data, SKNY performs two analyses: 1) a quantitative evaluation of the peripheral and intratumoral cells within a ROI, and 2) a "single tumor microenvironment analysis" with bulk gene expression within each contour.

AnnData object-based programming makes it compatible with scanpy, squidpy and stlearn.

* Free software: MIT licensed
* Documentation: https://skny.readthedocs.io.


.. image:: _images/SKYN_logo.svg
   :target: https://skny.readthedocs.io
   :width: 30%


Citation
--------

* TODO


