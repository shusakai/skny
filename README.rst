.. image:: https://img.shields.io/pypi/v/skny.svg
        :target: https://pypi.python.org/pypi/skny

.. image:: https://readthedocs.org/projects/skny/badge/?version=latest
        :target: https://skny.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

SKNY - SpatialKNifeY
=====================

.. image:: _images/SKYN_logo.svg
   :target: https://skny.readthedocs.io
   :width: 30%


**SKNY** is a tools for spatial anlysis stratified by distance from tumor for multiple platform such as `Xenium`_, `CosMx`_, and `PhenoCycler`_. 
It automatically contours the tumor based on the expression data and calculates the distance from these contours to each coordinate in space.
Using this distance data, SKNY performs two analyses: 

1. A quantitative evaluation of the peripheral and intratumoral cells.

2. "Single tumor microenvironment analysis" with bulk gene expression within each contour.

Also, `anndata`_ object-based programming makes it compatible with `scanpy`_, `squidpy`_ and `stlearn`_.


Tutorials
--------

* Documentation: https://skny.readthedocs.io.


Citation
--------

* TODO




.. _Xenium: https://www.10xgenomics.com/jp/platforms/xenium

.. _CosMx: https://nanostring.com/products/cosmx-spatial-molecular-imager/

.. _PhenoCycler: https://www.akoyabio.com/phenocycler/

.. _anndata: https://anndata.readthedocs.io/en/latest/

.. _scanpy: https://scanpy.readthedocs.io/en/stable/

.. _squidpy: https://squidpy.readthedocs.io/en/stable/

.. _stlearn: https://stlearn.readthedocs.io/en/latest/