"""Tests for skny.pp.calculate_region / region_density and the region plots (synthetic data)."""
import unittest

import numpy as np
import pandas as pd


def make_adata(seed=0):
    import anndata as ad
    rng = np.random.default_rng(seed)
    # tissue: 1500 x 1500 um; tumor nest = disc of radius 300 um centred at (700, 700)
    n = 20000
    xy = rng.uniform(0, 1500, size=(n, 2))
    r = np.hypot(xy[:, 0] - 700, xy[:, 1] - 700)
    label = np.where(r < 300, "Malignant", "Stroma").astype(object)
    immune = rng.random(n) < 0.1
    label[immune & (label == "Stroma")] = "T cell"
    # a few isolated mislabelled tumor cells far away (should not form a nest)
    label[(r > 600) & (rng.random(n) < 0.0005)] = "Malignant"
    obs = pd.DataFrame({"cell_type": pd.Categorical(label), "sample": "S1"}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(obs=obs)
    a.obsm["spatial"] = xy
    return a


class TestRegion(unittest.TestCase):
    def setUp(self):
        import skny as sk
        self.sk = sk
        self.a = make_adata()
        sk.pp.calculate_region(self.a, "Malignant", label_key="cell_type", sample_key="sample", verbose=False)

    def test_regions_assigned(self):
        a = self.a
        self.assertIn("region", a.obs)
        self.assertIn("region", a.uns)
        reg = a.obs["region"].astype(str)
        r = np.hypot(a.obsm["spatial"][:, 0] - 700, a.obsm["spatial"][:, 1] - 700)
        # cells deep inside the nest are core; cells 30-90 um outside are margin; cells far away are dist
        self.assertTrue((reg[r < 250] == "core").mean() > 0.99)
        self.assertTrue((reg[(r > 330) & (r < 390)] == "margin").mean() > 0.95)
        self.assertTrue((reg[r > 450] == "dist").mean() > 0.95)
        self.assertFalse((reg == "not_computed").any())

    def test_areas(self):
        area = self.a.uns["region"]["maps"]["S1"]["area_mm2"]
        self.assertAlmostEqual(area["core"], np.pi * 0.3 ** 2, delta=0.03)   # ~0.28 mm2
        self.assertGreater(area["dist"], 1.0)
        total = sum(area.values())            # tissue mask is smoothed, so it extends a little beyond the cells
        self.assertGreater(total, 1.5 * 1.5); self.assertLess(total, 1.5 * 1.5 * 1.5)

    def test_isolated_cells_not_core(self):
        a = self.a
        r = np.hypot(a.obsm["spatial"][:, 0] - 700, a.obsm["spatial"][:, 1] - 700)
        far_ref = (a.obs["cell_type"] == "Malignant").to_numpy() & (r > 600)
        self.assertGreater(far_ref.sum(), 0)
        self.assertTrue((a.obs.loc[far_ref, "region"].astype(str) == "dist").all())

    def test_margin_from_reference_legacy(self):
        a = make_adata()
        self.sk.pp.calculate_region(a, "Malignant", sample_key="sample", margin_from="reference", verbose=False)
        r = np.hypot(a.obsm["spatial"][:, 0] - 700, a.obsm["spatial"][:, 1] - 700)
        far_ref = (a.obs["cell_type"] == "Malignant").to_numpy() & (r > 600)
        # legacy definition: isolated reference cells sit in their own margin halo
        self.assertTrue((a.obs.loc[far_ref, "region"].astype(str) == "margin").all())
        self.assertGreater(a.uns["region"]["maps"]["S1"]["area_mm2"]["margin"], self.a.uns["region"]["maps"]["S1"]["area_mm2"]["margin"])

    def test_density(self):
        df = self.sk.pp.region_density(self.a, "cell_type", sample_key="sample")
        self.assertEqual(set(df["region"]), {"core", "margin", "dist"})
        t = df[(df["group"] == "T cell")].set_index("region")["density"]
        self.assertLess(t["core"], 0.2 * t["dist"])       # (almost) no T cells inside the nest
        # ~10% of 20000 cells on 2.25 mm2 (~890 cells/mm2), diluted by the smoothed tissue border
        self.assertGreater(t["dist"], 450); self.assertLess(t["dist"], 950)
        allc = df[df["group"] == "_all_cells"]
        self.assertEqual(int(allc["n"].sum()), int((self.a.obs["region"].astype(str).isin(["core", "margin", "dist"])).sum()))

    def test_skip_without_reference(self):
        a = make_adata(1)
        a.obs["cell_type"] = "Stroma"
        self.sk.pp.calculate_region(a, "Malignant", label_key="cell_type", verbose=False)
        self.assertTrue((a.obs["region"].astype(str) == "not_computed").all())
        self.assertEqual(len(a.uns["region"]["maps"]), 0)

    def test_coord_keys_and_copy(self):
        a = make_adata(2)
        a.obs["x_um"] = a.obsm["spatial"][:, 0]; a.obs["y_um"] = a.obsm["spatial"][:, 1]
        del a.obsm["spatial"]
        b = self.sk.pp.calculate_region(a, ["Malignant"], coord_keys=("x_um", "y_um"), copy=True, verbose=False)
        self.assertNotIn("region", a.obs)
        self.assertIn("region", b.obs)

    def test_plots(self):
        import matplotlib
        matplotlib.use("Agg")
        ax = self.sk.pl.region_map(self.a, sample="S1", highlight="T cell", show=False)
        self.assertIsNotNone(ax)
        df = self.sk.pp.region_density(self.a, "cell_type", sample_key="sample")
        fig, axes = self.sk.pl.region_density_plot(df, groups=["T cell"], min_area_mm2=0.2, show=False)
        self.assertEqual(len(axes), 3)


if __name__ == "__main__":
    unittest.main()
