"""SKNY (SpatialKNifeY): tumor microenvironment analysis for single-cell spatial omics.

Two workflows are provided:

* ``skny.pp.calculate_distance`` - expression-based tumor contour on a 10 um grid and the
  distance of every bin from the contour (the method of the paper), with
  ``skny.pl.distance_gene_plot``, ``skny.pl.clustering_heatmap``, ``skny.pl.deg_manhattan_plot``
  and per-tumor-nest pseudo-bulk (``skny.pp.convert_indivisual_solid`` / ``_peri``).
* ``skny.pp.calculate_region`` - tumor ``core`` / ``margin`` / ``dist`` regions from a cell
  label, ``skny.pp.region_density`` (cells / mm^2 per region), ``skny.pl.region_map`` and
  ``skny.pl.region_density_plot``.

Citation: Sakai SA et al., PLOS Comput Biol 21(2): e1012854 (2025),
https://doi.org/10.1371/journal.pcbi.1012854
"""

__author__ = """Shunsuke A. Sakai"""
__email__ = 'shusakai@east.ncc.go.jp'
__version__ = '0.3.2'
__url__  = 'https://github.com/shusakai/skny'


from . import pp
from . import pl
