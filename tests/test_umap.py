"""
Unit tests for Pure NumPy / SciPy UMAP implementation in PyChem-Pro.
Uses standard library unittest so no external pytest dependency is required.
"""
import unittest
import numpy as np

from src.algorithms.umap import UMAP
from src.algorithms.umap.metrics import pairwise_distances, find_knn
from src.algorithms.umap.fuzzy_simplicial_set import smooth_knn_dist, fuzzy_simplicial_set
from src.algorithms.umap.curve_fitting import find_ab_params
from src.features.chemical_space.analysis_engine import AnalysisEngine

class TestUMAP(unittest.TestCase):

    def setUp(self):
        np.random.seed(42)
        self.dummy_data = np.random.normal(size=(50, 10))
        self.binary_fingerprints = np.random.randint(0, 2, size=(40, 64)).astype(np.float64)

    def test_distance_metrics(self):
        # Test Euclidean distance
        D_euc = pairwise_distances(self.dummy_data, metric='euclidean')
        self.assertEqual(D_euc.shape, (50, 50))
        self.assertTrue(np.allclose(np.diag(D_euc), 0.0))

        # Test Tanimoto distance on binary fingerprints
        D_tan = pairwise_distances(self.binary_fingerprints, metric='tanimoto')
        self.assertEqual(D_tan.shape, (40, 40))
        self.assertTrue(np.all(D_tan >= 0.0) and np.all(D_tan <= 1.0))
        self.assertTrue(np.allclose(np.diag(D_tan), 0.0))

    def test_find_knn(self):
        D = pairwise_distances(self.dummy_data, metric='euclidean')
        knn_idx, knn_dists = find_knn(D, n_neighbors=10)
        self.assertEqual(knn_idx.shape, (50, 10))
        self.assertEqual(knn_dists.shape, (50, 10))
        self.assertTrue(np.all(knn_dists >= 0.0))

    def test_smooth_knn_dist(self):
        D = pairwise_distances(self.dummy_data, metric='euclidean')
        knn_idx, knn_dists = find_knn(D, n_neighbors=15)
        rho, sigma = smooth_knn_dist(knn_dists, k=15)
        self.assertEqual(rho.shape, (50,))
        self.assertEqual(sigma.shape, (50,))
        self.assertTrue(np.all(rho >= 0.0))
        self.assertTrue(np.all(sigma > 0.0))

    def test_curve_fitting(self):
        a, b = find_ab_params(spread=1.0, min_dist=0.1)
        self.assertIsInstance(a, float)
        self.assertIsInstance(b, float)
        self.assertGreater(a, 0.0)
        self.assertGreater(b, 0.0)

    def test_umap_fit_transform(self):
        umap_model = UMAP(n_neighbors=15, n_components=2, n_epochs=50, random_state=42)
        embedding = umap_model.fit_transform(self.dummy_data)
        
        self.assertEqual(embedding.shape, (50, 2))
        self.assertFalse(np.isnan(embedding).any())
        self.assertFalse(np.isinf(embedding).any())

    def test_umap_transform_out_of_sample(self):
        X_train = self.dummy_data[:40]
        X_test = self.dummy_data[40:]

        umap_model = UMAP(n_neighbors=10, n_components=2, n_epochs=30, random_state=42)
        umap_model.fit(X_train)

        Y_test = umap_model.transform(X_test)
        self.assertEqual(Y_test.shape, (10, 2))
        self.assertFalse(np.isnan(Y_test).any())

    def test_analysis_engine_integration(self):
        engine = AnalysisEngine()
        success, msg = engine.run_umap(self.dummy_data, n_neighbors=10, n_epochs=30)
        
        self.assertTrue(success)
        self.assertEqual(msg, "UMAP completed")
        self.assertIsNotNone(engine.embedding)
        self.assertEqual(engine.embedding.shape, (50, 2))

    def test_multiprocessing(self):
        umap_parallel = UMAP(n_neighbors=10, n_components=2, n_epochs=30, n_jobs=2, random_state=42)
        embedding = umap_parallel.fit_transform(self.dummy_data)
        
        self.assertEqual(embedding.shape, (50, 2))
        self.assertFalse(np.isnan(embedding).any())
        self.assertEqual(umap_parallel.n_jobs_resolved_, 2)

if __name__ == '__main__':
    unittest.main()
