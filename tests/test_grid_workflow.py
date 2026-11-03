"""End-to-end regression test of the expression-based grid workflow on synthetic data.

A square "tumor" expressing CDH1 is placed in a random cell cloud, the cells are binned into a
10 um grid (a minimal re-implementation of ``stlearn.tl.cci.grid`` so that stlearn is not
required), and the whole workflow 1 pipeline is run: calculate_distance -> plots ->
convert_indivisual_solid / convert_indivisual_peri.
"""
import unittest
import warnings

import numpy as np
import pandas as pd


def make_grid(seed=0, n=20000):
    import anndata as ad

    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 1000, n); y = rng.uniform(0, 800, n)
    tumor = (x > 300) & (x < 700) & (y > 250) & (y < 550)
    genes = ["CDH1", "KRT15", "ACTA2", "PTPRC", "COL1A1"]
    X = np.zeros((n, len(genes)))
    X[tumor, 0] = rng.poisson(5, tumor.sum()); X[tumor, 1] = rng.poisson(3, tumor.sum())
    X[~tumor, 2] = rng.poisson(4, (~tumor).sum()); X[~tumor, 4] = rng.poisson(3, (~tumor).sum())
    X[:, 3] = rng.poisson(1, n)
    n_col = int((x.max() - x.min()) / 10); n_row = int((y.max() - y.min()) / 10)
    xe = np.linspace(x.min(), x.max(), n_col + 1); ye = np.linspace(y.min(), y.max(), n_row + 1)
    cx = np.clip(np.digitize(x, xe) - 1, 0, n_col - 1); cy = np.clip(np.digitize(y, ye) - 1, 0, n_row - 1)
    df = pd.DataFrame(X, columns=genes); df["g"] = cx * n_row + cy
    agg = df.groupby("g").sum()
    grid = ad.AnnData(agg.values.astype(np.float32),
                      obs=pd.DataFrame(index=[f"grid_{i+1}" for i in agg.index]),
                      var=pd.DataFrame(index=genes))
    grid.uns["grid_xedges"] = xe; grid.uns["grid_yedges"] = ye
    return grid


class TestGridWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib
        matplotlib.use("Agg")
        import plotly.graph_objects as go
        go.Figure.show = lambda *a, **k: None  # never open a browser in tests
        warnings.simplefilter("ignore")
        import skny as sk
        cls.sk = sk
        cls.grid = sk.pp.calculate_distance(make_grid(), pos_marker_ls=["CDH1"])

    def test_calculate_distance(self):
        g = self.grid
        for key in ["marker", "marker_median", "marker_delineation", "shotest_30_delineation"]:
            self.assertIn(key, g.uns)
        sh = g.shortest
        self.assertIn("region", sh.columns)
        bands = sh.dropna()["region"]
        self.assertTrue((bands.apply(lambda i: i.right) < 0).any())   # inside the tumor
        self.assertTrue((bands.apply(lambda i: i.right) > 0).any())   # peri-tumoral stroma

    def test_plots(self):
        sk, g = self.sk, self.grid
        sk.pl.preprocessed_images(g, figsize=(600, None))
        sk.pl.equidistant_curves(g, figsize=(400, None))
        sk.pl.distance_gene_plot(g, "KRT15", show=False)
        sk.pl.distance_gene_plot(g, "KRT15", std=True, show=False)
        sk.pl.clustering_heatmap(g, std=True)
        sk.pl.deg_manhattan_plot(g, show=False)
        tbl = getattr(g, "gene_expresion_by_region")
        # KRT15 (tumor gene) must be higher inside the tumor than far outside
        self.assertGreater(tbl.loc["KRT15"].iloc[0], tbl.loc["KRT15"].iloc[-1])

    def test_single_tme(self):
        sk, g = self.sk, self.grid
        solid = sk.pp.convert_indivisual_solid(g)
        peri = sk.pp.convert_indivisual_peri(g, section=30)
        self.assertEqual(solid.n_obs, 1); self.assertEqual(peri.n_obs, 1)
        self.assertEqual(solid.var_names.tolist(), g.var_names.tolist())
        self.assertGreater(float(solid[:, "CDH1"].X.mean()), float(peri[:, "CDH1"].X.mean()))


if __name__ == "__main__":
    unittest.main()
