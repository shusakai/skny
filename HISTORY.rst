=======
History
=======

Unreleased
----------

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

