.. image:: https://img.shields.io/pypi/v/skny.svg
        :target: https://pypi.python.org/pypi/skny
        :alt: PyPI

.. image:: https://readthedocs.org/projects/skny/badge/?version=latest
        :target: https://skny.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/badge/paper-PLOS%20Comput%20Biol%202025-blue.svg
        :target: https://doi.org/10.1371/journal.pcbi.1012854
        :alt: Paper

.. image:: https://img.shields.io/github/license/shusakai/skny.svg
        :target: https://github.com/shusakai/skny/blob/main/LICENSE
        :alt: License

SKNY (SpatialKNifeY) - tumor microenvironment analysis for single-cell spatial omics
====================================================================================

.. image:: _images/SKYN_logo.svg
   :target: https://skny.readthedocs.io
   :width: 30%

**SKNY** (SpatialKNifeY) is a Python package for **spatial transcriptomics / spatial proteomics
analysis of the tumor microenvironment (TME)**. It works on single-cell spatial omics data from
`Xenium`_, `CosMx`_, `PhenoCycler`_ and any other platform that can be loaded into an `anndata`_
object, so it fits into `scanpy`_, `squidpy`_ and `stlearn`_ workflows.

SKNY answers questions such as:

* Where is the **tumor border** in this section, and how far is each location or cell from it?
* How does gene expression change with **distance from the tumor** (intra-tumoral vs peri-tumoral)?
* Which genes are **differentially expressed between the tumor core and the surrounding stroma**?
* What is the **density of immune cells (cells / mm²) in the tumor core, the tumor margin and distant tissue**?
* How do individual tumor nests ("single tumor microenvironment") differ from each other?

Published in *PLOS Computational Biology* (2025): https://doi.org/10.1371/journal.pcbi.1012854
(see `Citation`_).


Two workflows
-------------

**1. Expression-based tumor contour on a grid** (``sk.pp.calculate_distance``; the method in the paper)

The tissue is binned into a 10 µm grid, the tumor is contoured from marker gene expression
(e.g. ``CDH1``, ``EPCAM``, keratins) with image processing, and the shortest-path distance from
the contour to every grid bin is calculated. Downstream functions then give:

* gene expression as a function of distance from the tumor border (``sk.pl.distance_gene_plot``),
* clustering of all genes by their distance profile (``sk.pl.clustering_heatmap``),
* differentially expressed genes between distance bands (``sk.pl.deg_manhattan_plot``),
* per-tumor-nest ("single-TME") pseudo-bulk expression inside each nest
  (``sk.pp.convert_indivisual_solid``) and in the band around it (``sk.pp.convert_indivisual_peri``).

**2. Cell-label based regions: tumor core / margin / dist** (``sk.pp.calculate_region``; since v0.3.0)

For cell-level data with annotated cell types (scanpy / Xenium 5K / CosMx / PhenoCycler), the tumor
is defined by a cell label such as ``"Malignant"``. Every cell is assigned to the tumor ``core``
(within 15 µm of a tumor cell), the peri-tumoral ``margin`` (stroma within 100 µm of the core) or
``dist`` (farther tissue), using an all-cell tissue mask so that empty space is never counted.
``sk.pp.region_density`` then reports the density (cells / mm²) of any cell group per region, and
``sk.pl.region_map`` / ``sk.pl.region_density_plot`` draw the results.


Installation
------------

.. code-block:: console

   pip install skny

SKNY depends on ``stlearn`` (for the grid workflow), ``scanpy``/``anndata``, ``opencv-python``,
``networkx``, ``plotly`` and ``matplotlib``. We currently recommend a dedicated environment with
Python 3.8 to 3.10 because of the pinned ``numpy`` / ``numba`` versions required by ``stlearn``:

.. code-block:: console

   conda create -n skny python=3.10
   conda activate skny
   pip install skny


Quick start
-----------

**Workflow 1: distance from the expression-based tumor contour (Xenium example)**

.. code-block:: python

   import skny as sk
   import stlearn as st

   # Load a Xenium run with stlearn (any AnnData with obs["imagecol"], obs["imagerow"] works)
   adata = st.ReadXenium(
       feature_cell_matrix_file="Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5",
       cell_summary_file="Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz",
       library_id="example", image_path=None, scale=1, spot_diameter_fullres=10,
   )

   # Bin cells into a 10 um grid
   n_col = int((adata.obs.imagecol.max() - adata.obs.imagecol.min()) / 10)
   n_row = int((adata.obs.imagerow.max() - adata.obs.imagerow.min()) / 10)
   grid = st.tl.cci.grid(adata, n_row=n_row, n_col=n_col, n_cpus=8, verbose=False)

   # Contour the tumor from marker genes and compute the distance of every bin from the contour
   grid = sk.pp.calculate_distance(grid, pos_marker_ls=["CDH1"])

   sk.pl.preprocessed_images(grid, figsize=(800, None))   # tumor mask, contour, equidistant curves
   sk.pl.distance_gene_plot(grid, "KRT15")                 # expression vs distance from the tumor
   sk.pl.clustering_heatmap(grid, std=True)                # all genes clustered by distance profile
   sk.pl.deg_manhattan_plot(grid)                          # DEGs between intra- and peri-tumoral bands

   # "Single tumor microenvironment": pseudo-bulk expression per tumor nest
   solid = sk.pp.convert_indivisual_solid(grid)            # inside each nest  (-inf, 0] um
   peri  = sk.pp.convert_indivisual_peri(grid, section=30) # 0 to 30 um band around each nest

**Workflow 2: tumor core / margin / dist from a cell label**

.. code-block:: python

   import skny as sk

   # adata: one row per cell, coordinates in micrometres in adata.obsm["spatial"],
   #        cell types in adata.obs["cell_type"], sample id in adata.obs["sample"]
   sk.pp.calculate_region(adata, "Malignant", label_key="cell_type", sample_key="sample")
   # -> adata.obs["region"] in {"core", "margin", "dist", "outside"}, adata.obs["region_dist_um"]

   df = sk.pp.region_density(adata, "immune_subtype", sample_key="sample")
   # -> columns: sample, region, group, n, area_mm2, density (cells / mm^2)

   sk.pl.region_map(adata, sample="S1", highlight=["T cell", "Myeloid"])
   sk.pl.region_density_plot(df, groups=["CD8+ T", "Treg", "M2-like"])


When to use SKNY (and how it relates to scanpy / squidpy)
---------------------------------------------------------

* `scanpy`_ handles clustering and differential expression; `squidpy`_ handles neighborhood graphs,
  co-occurrence and ligand-receptor analysis. Neither defines a **tumor border** or measures
  **distance from it**.
* Use SKNY when your question is about the **spatial relationship to the tumor**: expression
  gradients from the tumor border, immune infiltration into the tumor core vs the margin,
  comparison of individual tumor nests, or anything stratified by distance from a tumor domain.
* SKNY takes and returns ``AnnData`` objects, so you can keep using scanpy / squidpy before and
  after it (e.g. run ``sc.tl.rank_genes_groups`` on ``adata.obs["region"]``).


Documentation and tutorials
---------------------------

* Documentation: https://skny.readthedocs.io
* Tutorials (Jupyter notebooks, public 10x Xenium breast cancer data):

  * `Xenium analysis`_ - workflow 1 end to end
  * `Single-TME analysis`_ - per-tumor-nest analysis
  * `Region analysis`_ - workflow 2 (core / margin / dist, cell density)
  * `Reproduction of the paper figures`_

* API reference: https://skny.readthedocs.io/en/latest/api.html
* Machine-readable summary for coding assistants: `llms.txt`_ and `AGENTS.md`_


Citation
--------

If you use SKNY, please cite:

  Sakai SA, Nomura R, Nagasawa S, Chi S, Suzuki A, Suzuki Y, Imai M, Nakamura Y, Yoshino T,
  Ishikawa S, Tsuchihara K, Kageyama S-I, Yamashita R.
  **SpatialKNifeY (SKNY): Extending from spatial domain to surrounding area to identify
  microenvironment features with single-cell spatial omics data.**
  *PLOS Computational Biology* 21(2): e1012854 (2025).
  https://doi.org/10.1371/journal.pcbi.1012854

.. code-block:: bibtex

   @article{Sakai2025SKNY,
     title   = {SpatialKNifeY (SKNY): Extending from spatial domain to surrounding area to identify microenvironment features with single-cell spatial omics data},
     author  = {Sakai, Shunsuke A. and Nomura, Ryosuke and Nagasawa, Satoi and Chi, SungGi and Suzuki, Ayako and Suzuki, Yutaka and Imai, Mitsuho and Nakamura, Yoshiaki and Yoshino, Takayuki and Ishikawa, Shumpei and Tsuchihara, Katsuya and Kageyama, Shun-Ichiro and Yamashita, Riu},
     journal = {PLOS Computational Biology},
     year    = {2025},
     volume  = {21},
     number  = {2},
     pages   = {e1012854},
     doi     = {10.1371/journal.pcbi.1012854}
   }

A ``CITATION.cff`` file is included in the repository (GitHub "Cite this repository" button).


License
-------

MIT. See `LICENSE`_.


.. _Xenium: https://www.10xgenomics.com/platforms/xenium
.. _CosMx: https://nanostring.com/products/cosmx-spatial-molecular-imager/
.. _PhenoCycler: https://www.akoyabio.com/phenocycler/
.. _anndata: https://anndata.readthedocs.io/en/latest/
.. _scanpy: https://scanpy.readthedocs.io/en/stable/
.. _squidpy: https://squidpy.readthedocs.io/en/stable/
.. _stlearn: https://stlearn.readthedocs.io/en/latest/
.. _Xenium analysis: https://skny.readthedocs.io/en/latest/notebooks/Xenium_analysis.html
.. _Single-TME analysis: https://skny.readthedocs.io/en/latest/notebooks/single-TME_analysis.html
.. _Region analysis: https://skny.readthedocs.io/en/latest/notebooks/region_analysis.html
.. _Reproduction of the paper figures: https://skny.readthedocs.io/en/latest/notebooks/SKNY_paper_Sakai_et_al.html
.. _llms.txt: https://github.com/shusakai/skny/blob/main/llms.txt
.. _AGENTS.md: https://github.com/shusakai/skny/blob/main/AGENTS.md
.. _LICENSE: https://github.com/shusakai/skny/blob/main/LICENSE
