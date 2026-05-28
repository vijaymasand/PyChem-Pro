"""
Chemical Space Visualization - Analysis Engine
Handles dimensionality reduction and clustering using scikit-learn.
"""
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN
from sklearn.ensemble import IsolationForest
import numpy as np
import logging

class AnalysisEngine:
    def __init__(self):
        self.embedding = None
        self.clusters = None
        self.outliers = None

    def run_pca(self, X: np.ndarray, n_components: int = 2):
        """Runs PCA embedding."""
        try:
            pca = PCA(n_components=n_components)
            self.embedding = pca.fit_transform(X)
            return True, "PCA completed"
        except Exception as e:
            logging.error(f"PCA Error: {e}")
            return False, str(e)

    def run_tsne(self, X: np.ndarray, perplexity: float = 30.0, learning_rate: str = 'auto', random_state: int = 42):
        """Runs t-SNE embedding."""
        try:
            # handle 'auto' string conversion or keep it if sklearn supports it
            lr = learning_rate if isinstance(learning_rate, str) else float(learning_rate)
            tsne = TSNE(n_components=2, perplexity=perplexity, learning_rate=lr, random_state=random_state, init='pca')
            self.embedding = tsne.fit_transform(X)
            return True, "t-SNE completed"
        except Exception as e:
            logging.error(f"t-SNE Error: {e}")
            return False, str(e)

    def run_kmeans(self, n_clusters: int = 5):
        """Runs K-Means clustering on the embeddings."""
        if self.embedding is None:
            return False, "No embeddings available for clustering"
        try:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            self.clusters = kmeans.fit_predict(self.embedding)
            return True, "K-Means completed"
        except Exception as e:
            logging.error(f"K-Means Error: {e}")
            return False, str(e)

    def run_dbscan(self, eps: float = 0.5, min_samples: int = 5):
        """Runs DBSCAN clustering on the embeddings."""
        if self.embedding is None:
            return False, "No embeddings available for clustering"
        try:
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            self.clusters = dbscan.fit_predict(self.embedding)
            return True, "DBSCAN completed"
        except Exception as e:
            logging.error(f"DBSCAN Error: {e}")
            return False, str(e)

    def detect_outliers(self, X: np.ndarray, contamination: float = 0.05):
        """Runs Isolation Forest on the original features X."""
        try:
            iso = IsolationForest(contamination=contamination, random_state=42)
            # -1 for outliers, 1 for inliers
            preds = iso.fit_predict(X)
            self.outliers = preds == -1
            return True, "Outlier detection completed"
        except Exception as e:
            logging.error(f"Isolation Forest Error: {e}")
            return False, str(e)
