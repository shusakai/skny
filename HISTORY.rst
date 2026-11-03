=======
History
=======

0.3.4 (2026-11-03)
------------------

* Configurable image processing in ``skny.pp.calculate_distance``: ``filters`` (a list of steps:
  median / box / gaussian filters, morphological open / close / erode / dilate, hole filling,
  minimum component area, explicit threshold) and ``threshold``. The default
  ``[("median", 3), ("box", 3)]`` with threshold 90 is the setting of the paper. The settings are
  recorded in ``uns["image_processing"]``; ``skny.pp.apply_filters`` previews a filter set.
* The smoothed mask now defines both the contour and the tumor pixels (previously the box
  filter + threshold was used for the contour only, the median-filtered image for the pixels).
  On the 10x Xenium breast cancer section this changes the intra-tumoral bins from 66.8k to
  69.5k (+4 %) and the number of tumor nests from 431 to 535 with the default setting.
* At least one smoothing step is now required (without smoothing the contour fragments into
  thousands of pieces and the distance calculation takes tens of minutes).

0.3.3 (2026-11-03)
------------------

* ``skny.pp.calculate_distance`` can now define the tumor spatial domain from **annotated cell
  labels** (``use_label`` + ``pos_label_ls`` / ``neg_label_ls``, with ``min_cells`` and
  ``min_fraction``), in addition to marker genes. Build the grid with
  ``skny.pp.grid(adata, use_label="cell_type")``. The definition used is stored in
  ``uns["tumor_definition"]``.
* **Tumor interior fix.** The binary tumor image is now built directly from the bins instead of
  being rendered through a matplotlib figure. The rendering interpolated the colours, so most
  tumor bins no longer had the exact colour the pipeline looks for and were treated as *outside*
  the tumor: on the 10x Xenium breast cancer section only 13k of 78k CDH1+ bins survived
  (3.5k after the median filter) and 27k bins got a negative (intra-tumoral) distance, versus
  67k now. Intra-tumoral bands and per-nest ("solid") objects are therefore larger than with
  <= 0.3.2; the contour itself and the peri-tumoral bands change little.

0.3.2 (2026-11-03)
------------------

* New ``skny.pp.read_xenium`` (10x Xenium reader, h5py only) and ``skny.pp.grid`` (10 um
  gridding with the same output format as ``stlearn.tl.cci.grid``), so the expression-based
  workflow runs with ``pip install skny`` alone; ``skny.pl.gene_plot`` replaces
  ``stlearn.pl.gene_plot`` in the tutorials.
* Fixed a one-bin offset: stlearn numbers grid bins from ``grid_0`` while SKNY assumed
  ``grid_1``, which shifted every mask by one bin (10 um). The numbering is now detected
  (``skny.pp.bin_names``); results on stlearn grids move by one bin compared with <= 0.3.1.
* Tutorials no longer depend on stlearn (the *Tumor contour with Xenium data* notebook was
  re-executed with skny 0.3.2); the sdist no longer ships the notebooks (11 MB -> < 1 MB).

0.3.1 (2026-11-03)
------------------

* Dependencies relaxed: no more pins on ``numpy`` / ``numba``; ``stlearn`` and ``scanpy`` are no
  longer required (install ``skny[stlearn]`` for the grid tutorial). SKNY now runs on
  Python 3.8 to 3.13 with numpy 2 and pandas 3.
* Compatibility fixes for recent matplotlib / pandas (figure buffer access, row-wise z-score,
  ``applymap`` removal, NaN-safe distance-band handling in ``convert_indivisual_*``).
* New end-to-end regression test of the grid workflow on synthetic data
  (``tests/test_grid_workflow.py``, no external data needed).

* Citation information (PLOS Computational Biology 2025) added to the README and a
  ``CITATION.cff`` file; ``llms.txt`` and ``AGENTS.md`` added for coding assistants.
* PyPI metadata (keywords, project URLs, classifiers) and NumPy-style docstrings for all
  public functions.
* Removed the unused ``skny`` console script (template leftover).
* ``skny.pp.calculate_region`` now warns instead of raising when no cell carries the
  reference label (all cells become ``not_computed``), as covered by the unit tests.

0.3.0 (2026-11-02)
------------------

* New cell-label based region analysis (a variant of ``calculate_distance`` that works on
  cell-level AnnData instead of marker expression on a grid):
  ``skny.pp.calculate_region`` (tumor ``core`` / ``margin`` / ``dist`` from the Euclidean
  distance to the nearest reference cell, with an all-cell tissue mask),
  ``skny.pp.region_density`` (cells / mm^2 of any cell group per region),
  ``skny.pl.region_map`` and ``skny.pl.region_density_plot``.
* Tutorial ``notebooks/region_analysis`` and unit tests ``tests/test_region.py``.


0.2.16 (2024-02-14)
-----------------------
0.2.8 (2024-02-01)
------------------
0.2.7 (2024-01-31)
------------------
0.2.6 (2024-01-31)
------------------

* Practice of upload to PyPI.

0.2.5 (2024-01-30)
------------------
0.2.4 (2024-01-30)
------------------
0.2.3 (2024-01-30)
------------------
0.2.2 (2024-01-30)
------------------
0.2.1 (2024-01-30)
------------------
0.1.2 (2024-01-30)
------------------
0.1.1 (2024-01-30)
------------------
0.1.0 (2024-01-21)
------------------

* First release on PyPI.

