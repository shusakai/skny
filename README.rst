.. image:: https://img.shields.io/pypi/v/skny.svg
        :target: https://pypi.python.org/pypi/skny

.. image:: https://readthedocs.org/projects/skny/badge/?version=latest
        :target: https://skny.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

SKNY - Spatial omics analysis tools for tumor microenvironment 
==================================================================================== 

.. image:: _images/SKYN_logo.svg
   :target: https://skny.readthedocs.io
   :width: 30%


**SKNY** (SpatialKNifeY) is a tools for spatial omics analysis of tumor microenvironment for multiple platform such as `Xenium`_, `CosMx`_, and `PhenoCycler`_. 
SKNY automatically contours the tumor based on the expression data and calculates the distance from these contours to each coordinate in space.
Using this distance data, we provide two analyses: 

1. A quantitative evaluation of the peri- and intra-tumoral gene expressions.

2. "Single tumor microenvironment analysis" with bulk gene expression within each contour.

Since v0.3.0 SKNY also provides a **cell-label based variant** (``skny.pp.calculate_region``) for
cell-level data with annotated cell types (e.g. scanpy / Xenium 5K): the tumor is defined by a cell
label such as ``"Malignant"`` and every cell is assigned to the tumor ``core`` (within 15 um of a
tumor cell), the peri-tumoral ``margin`` (stroma within 100 um of the core) or ``dist`` (farther
tissue), with an all-cell tissue mask so that empty space is never counted. ``skny.pp.region_density``
then gives the density (cells / mm^2) of any cell group (e.g. immune subtypes) per region, and
``skny.pl.region_map`` / ``skny.pl.region_density_plot`` draw the results.

.. code-block:: python

   import skny as sk
   sk.pp.calculate_region(adata, "Malignant", label_key="cell_type", sample_key="sample")  # coordinates in um in adata.obsm["spatial"]
   df = sk.pp.region_density(adata, "immune_subtype", sample_key="sample")
   sk.pl.region_map(adata, sample="S1", highlight=["T cell", "Myeloid"])
   sk.pl.region_density_plot(df, groups=["CD8+ T", "Treg", "M2-like"])

Also, `anndata`_ object-based programming makes it compatible with `scanpy`_, `squidpy`_ and `stlearn`_.

**Recommended Python version:**  
We recommend using Python 3.8 or 3.9 to ensure compatibility.


Documentation
----------------

* https://skny.readthedocs.io.


Citation
----------------

* TODO




.. _Xenium: https://www.10xgenomics.com/jp/platforms/xenium

.. _CosMx: https://nanostring.com/products/cosmx-spatial-molecular-imager/

.. _PhenoCycler: https://www.akoyabio.com/phenocycler/

.. _anndata: https://anndata.readthedocs.io/en/latest/

.. _scanpy: https://scanpy.readthedocs.io/en/stable/

.. _squidpy: https://squidpy.readthedocs.io/en/stable/

.. _stlearn: https://stlearn.readthedocs.io/en/latest/
