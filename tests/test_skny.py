#!/usr/bin/env python

"""Smoke tests for the `skny` package: import and public API."""

import unittest


class TestSknyPackage(unittest.TestCase):
    def test_import_and_version(self):
        import skny

        self.assertIsInstance(skny.__version__, str)

    def test_public_api(self):
        import skny as sk

        for name in [
            "calculate_distance",
            "convert_indivisual_solid",
            "convert_indivisual_peri",
            "calculate_region",
            "regions_from_coords",
            "region_density",
        ]:
            self.assertTrue(callable(getattr(sk.pp, name)), name)

        for name in [
            "preprocessed_images",
            "equidistant_curves",
            "clustering_heatmap",
            "distance_gene_plot",
            "deg_manhattan_plot",
            "region_map",
            "region_density_plot",
        ]:
            self.assertTrue(callable(getattr(sk.pl, name)), name)


if __name__ == "__main__":
    unittest.main()
