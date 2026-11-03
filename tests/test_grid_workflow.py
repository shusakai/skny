"""End-to-end regression tests of the expression-based grid workflow on synthetic data.

A square "tumor" expressing CDH1 is placed in a random cell cloud, the cells are binned into a
10 um grid with ``skny.pp.grid`` (no stlearn needed), and the whole workflow 1 pipeline is run:
calculate_distance -> plots -> convert_indivisual_solid / convert_indivisual_peri. The
``skny.pp.read_xenium`` reader is tested on a small synthetic 10x HDF5 file.
"""
import unittest
import warnings

import numpy as np
import pandas as pd


def make_cells(seed=0, n=20000):
    """Random cell cloud (1000 x 800 um) with a square tumor expressing CDH1 / KRT15."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 1000, n); y = rng.uniform(0, 800, n)
    tumor = (x > 300) & (x < 700) & (y > 250) & (y < 550)
    genes = ["CDH1", "KRT15", "ACTA2", "PTPRC", "COL1A1"]
    X = np.zeros((n, len(genes)))
    X[tumor, 0] = rng.poisson(5, tumor.sum()); X[tumor, 1] = rng.poisson(3, tumor.sum())
    X[~tumor, 2] = rng.poisson(4, (~tumor).sum()); X[~tumor, 4] = rng.poisson(3, (~tumor).sum())
    X[:, 3] = rng.poisson(1, n)
    adata = ad.AnnData(X.astype(np.float32),
                       obs=pd.DataFrame({"imagecol": x, "imagerow": y, "cell_type": np.where(tumor, "Tumor", "Stroma")},
                                        index=[f"c{i}" for i in range(n)]),
                       var=pd.DataFrame(index=genes))
    adata.obsm["spatial"] = np.c_[x, y]
    return adata


def make_grid(seed=0, n=20000, one_based=False):
    import skny as sk

    grid = sk.pp.grid(make_cells(seed, n), bin_um=10, verbose=False)
    if one_based:  # emulate a 1-based object (no offset marker; detected from the bin centres)
        grid.obs_names = [f"grid_{int(n.split('_')[1]) + 1}" for n in grid.obs_names]
        del grid.uns["grid_bin_offset"]
    return grid


class TestGridAndReader(unittest.TestCase):
    def test_grid_matches_stlearn_convention(self):
        import skny as sk

        cells = make_cells(n=5000)
        g = sk.pp.grid(cells, bin_um=10, use_label="cell_type", verbose=False)
        n_row, n_col = sk.pp.grid_shape(g)
        self.assertEqual((n_row, n_col), (len(g.uns["grid_yedges"]) - 1, len(g.uns["grid_xedges"]) - 1))
        self.assertTrue(g.obs_names[0].startswith("grid_"))
        self.assertEqual(g.uns["grid_counts"].shape, (n_col, n_row))
        self.assertEqual(int(g.obs["n_cells"].sum()), cells.n_obs)
        # bin n = col * n_row + row and the bin centre must contain the cells summed into it
        nums = np.array([int(s.split("_")[1]) for s in g.obs_names])
        col, row = nums // n_row, nums % n_row
        np.testing.assert_allclose(g.obs["imagecol"].to_numpy(), (g.uns["grid_xedges"][col] + g.uns["grid_xedges"][col + 1]) / 2)
        np.testing.assert_allclose(g.obs["imagerow"].to_numpy(), (g.uns["grid_yedges"][row] + g.uns["grid_yedges"][row + 1]) / 2)
        np.testing.assert_allclose(np.asarray(g.X).sum(), np.asarray(cells.X).sum(), rtol=1e-5)
        self.assertIn("cell_type", g.uns); self.assertIn("cell_type", g.obs)
        # label fractions sum to 1 for every non-empty bin
        np.testing.assert_allclose(g.uns["cell_type"].sum(1).to_numpy(), 1.0)
        # explicit n_row / n_col like stlearn
        g2 = sk.pp.grid(cells, n_row=40, n_col=50, verbose=False)
        self.assertEqual(sk.pp.grid_shape(g2), (40, 50))
        # obsm-based coordinates
        g3 = sk.pp.grid(cells, bin_um=10, spatial_key="spatial", verbose=False)
        self.assertEqual(g3.shape, g.shape)

    def test_label_based_contour(self):
        import skny as sk

        cells = make_cells(n=20000)
        g_lab = sk.pp.grid(cells, bin_um=10, use_label="cell_type", verbose=False)
        g_mrk = sk.pp.grid(cells, bin_um=10, verbose=False)
        d_lab = sk.pp.calculate_distance(g_lab, use_label="cell_type", pos_label_ls=["Tumor"])
        d_mrk = sk.pp.calculate_distance(g_mrk, pos_marker_ls=["CDH1"])
        self.assertEqual(d_lab.uns["tumor_definition"]["mode"], "label")
        self.assertEqual(d_mrk.uns["tumor_definition"]["mode"], "marker")
        # CDH1 is expressed only by "Tumor" cells, so both definitions give (almost) the same mask
        m1 = d_lab.uns["marker"][..., 0] > 0; m2 = d_mrk.uns["marker"][..., 0] > 0
        agree = (m1 == m2).mean()
        self.assertGreater(agree, 0.98)
        bands = d_lab.shortest.dropna()["region"].apply(lambda i: i.right)
        self.assertTrue((bands < 0).any() and (bands > 0).any())
        # stricter domain: min_fraction / min_cells shrink the mask, neg labels remove bins
        d_strict = sk.pp.calculate_distance(sk.pp.grid(cells, bin_um=10, use_label="cell_type", verbose=False),
                                            use_label="cell_type", pos_label_ls="Tumor", min_cells=2, min_fraction=0.5)
        self.assertLessEqual((d_strict.uns["marker"][..., 0] > 0).sum(), m1.sum())
        d_neg = sk.pp.calculate_distance(sk.pp.grid(cells, bin_um=10, use_label="cell_type", verbose=False),
                                         use_label="cell_type", pos_label_ls=["Tumor"], neg_label_ls=["Stroma"])
        self.assertLess((d_neg.uns["marker"][..., 0] > 0).sum(), m1.sum())
        # the downstream functions work on a label-based grid
        sk.pl.distance_gene_plot(d_lab, "KRT15", show=False)
        solid = sk.pp.convert_indivisual_solid(d_lab); self.assertEqual(solid.n_obs, 1)
        # errors: both / neither definitions, missing use_label, unknown label
        with self.assertRaises(ValueError):
            sk.pp.calculate_distance(g_mrk, pos_marker_ls=["CDH1"], pos_label_ls=["Tumor"])
        with self.assertRaises(ValueError):
            sk.pp.calculate_distance(g_mrk)
        with self.assertRaises(ValueError):
            sk.pp.calculate_distance(g_lab, pos_label_ls=["Tumor"])
        with self.assertRaises(KeyError):
            sk.pp.calculate_distance(g_mrk, use_label="cell_type", pos_label_ls=["Tumor"])
        with self.assertRaises(KeyError):
            sk.pp.calculate_distance(g_lab, use_label="cell_type", pos_label_ls=["Malignant"])

    def test_filters(self):
        import skny as sk

        cells = make_cells(n=20000)
        g = sk.pp.grid(cells, bin_um=10, verbose=False)
        raw = (g.to_df()["CDH1"] >= 1)
        # apply_filters on a boolean / colour image; default = paper setting
        d_def = sk.pp.calculate_distance(sk.pp.grid(cells, bin_um=10, verbose=False), pos_marker_ls=["CDH1"])
        self.assertEqual(d_def.uns["image_processing"], {"filters": ["median 3", "box 3"], "threshold": 90})
        # the processed grid can be written to / read from h5ad (uns must be serialisable)
        import tempfile, os, anndata as ad
        with tempfile.TemporaryDirectory() as tmp:
            d_def.write_h5ad(os.path.join(tmp, "g.h5ad")); back = ad.read_h5ad(os.path.join(tmp, "g.h5ad"))
        self.assertEqual(list(back.uns["image_processing"]["filters"]), ["median 3", "box 3"])
        self.assertEqual(back.uns["tumor_definition"]["mode"], "marker")
        m_def = sk.pp.apply_filters(d_def.uns["marker"])
        np.testing.assert_array_equal(m_def, d_def.uns["marker_median"][..., 0] > 0)
        # no smoothing keeps the raw mask
        m_none = sk.pp.apply_filters(d_def.uns["marker"], filters=None)
        np.testing.assert_array_equal(m_none, d_def.uns["marker"][..., 0] > 0)
        # alternative pipelines run end to end and give a sensible domain
        for filters, thr in [([("gaussian", 5), ("close", 3), "fill_holes"], 128),
                             ([("median", 5), ("min_area", 10)], 90),
                             (["median", "box", "threshold", ("open", 3)], 90),
                             ([{"op": "gaussian", "ksize": 3, "sigma": 1.0}], 90)]:
            d = sk.pp.calculate_distance(sk.pp.grid(cells, bin_um=10, verbose=False), pos_marker_ls=["CDH1"],
                                         filters=filters, threshold=thr)
            m = d.uns["marker_median"][..., 0] > 0
            self.assertGreater(m.sum(), 0.5 * m_def.sum()); self.assertLess(m.sum(), 1.5 * m_def.sum())
            self.assertTrue((d.shortest.dropna()["region"].apply(lambda i: i.right) < 0).any())
        # min_area removes an isolated small nest
        mask = np.zeros((40, 60), bool); mask[10:30, 10:40] = True; mask[2:4, 50:52] = True
        self.assertTrue(sk.pp.apply_filters(mask, filters=None)[2, 50])
        self.assertFalse(sk.pp.apply_filters(mask, filters=[("min_area", 10)])[2, 50])
        self.assertTrue(sk.pp.apply_filters(mask, filters=[("min_area", 10)])[20, 20])
        # fill_holes fills an enclosed hole
        mask2 = mask.copy(); mask2[18:22, 20:24] = False
        self.assertFalse(sk.pp.apply_filters(mask2, filters=None)[20, 22])
        self.assertTrue(sk.pp.apply_filters(mask2, filters=["fill_holes"])[20, 22])
        # errors
        with self.assertRaises(ValueError):
            sk.pp.calculate_distance(sk.pp.grid(cells, bin_um=10, verbose=False), pos_marker_ls=["CDH1"], filters=None)
        with self.assertRaises(ValueError):
            sk.pp.apply_filters(mask, filters=[("sharpen", 3)])
        with self.assertRaises(ValueError):
            sk.pp.apply_filters(mask, filters=[("median", 4)])

    def test_bin_names_zero_and_one_based(self):
        import skny as sk

        g0 = make_grid(n=3000); g1 = make_grid(n=3000, one_based=True)
        n_row, n_col = sk.pp.grid_shape(g0)
        self.assertEqual(sk.pp.bin_names(g0)[0], "grid_0"); self.assertEqual(len(sk.pp.bin_names(g0)), n_row * n_col)
        self.assertEqual(sk.pp.bin_names(g1)[0], "grid_1"); self.assertEqual(sk.pp.bin_names(g1)[-1], f"grid_{n_row * n_col}")
        # both conventions must give the same tumor mask after calculate_distance
        d0 = sk.pp.calculate_distance(g0, pos_marker_ls=["CDH1"]); d1 = sk.pp.calculate_distance(g1, pos_marker_ls=["CDH1"])
        np.testing.assert_array_equal(d0.uns["marker"], d1.uns["marker"])

    def test_read_xenium(self):
        import tempfile, os, h5py, gzip
        import scipy.sparse as sp
        import skny as sk

        cells = make_cells(n=500)
        X = sp.csc_matrix(np.asarray(cells.X).T)  # 10x stores features x barcodes in CSC
        with tempfile.TemporaryDirectory() as d:
            h5 = os.path.join(d, "test_cell_feature_matrix.h5")
            with h5py.File(h5, "w") as f:
                m = f.create_group("matrix")
                m["data"] = X.data.astype(np.int32); m["indices"] = X.indices.astype(np.int64); m["indptr"] = X.indptr.astype(np.int64)
                m["shape"] = np.array(X.shape, dtype=np.int32)  # (features, barcodes)
                m["barcodes"] = np.array([str(i + 1) for i in range(cells.n_obs)], dtype="S")
                feat = m.create_group("features")
                names = list(cells.var_names) ; names[-1] = "NegControl"
                feat["name"] = np.array(names, dtype="S"); feat["id"] = np.array([f"ENSG{i}" for i in range(len(names))], dtype="S")
                feat["feature_type"] = np.array(["Gene Expression"] * (len(names) - 1) + ["Negative Control Probe"], dtype="S")
            tbl = pd.DataFrame({"cell_id": [str(i + 1) for i in range(cells.n_obs)],
                                "x_centroid": cells.obs["imagecol"].to_numpy(), "y_centroid": cells.obs["imagerow"].to_numpy(),
                                "transcript_counts": np.asarray(cells.X).sum(1), "cell_area": 50.0})
            tbl.to_csv(os.path.join(d, "test_cells.csv.gz"), index=False)
            a = sk.pp.read_xenium(d, verbose=False)                        # directory form
            b = sk.pp.read_xenium(h5, cells_file=os.path.join(d, "test_cells.csv.gz"), genes_only=False, verbose=False)
        self.assertEqual(a.n_vars, cells.n_vars - 1); self.assertEqual(b.n_vars, cells.n_vars)
        self.assertEqual(a.n_obs, cells.n_obs)
        np.testing.assert_allclose(a.obs["imagecol"].to_numpy(), cells.obs["imagecol"].to_numpy())
        np.testing.assert_allclose(a[:, "CDH1"].X.toarray().ravel(), np.asarray(cells[:, "CDH1"].X).ravel())
        self.assertIn("spatial", a.obsm); self.assertIn("cell_area", a.obs)


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
        sk.pl.gene_plot(g, "KRT15"); sk.pl.gene_plot(g, ["KRT15", "CDH1"], zoom_coord=[200, 800, 200, 600], vmax=10)
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
