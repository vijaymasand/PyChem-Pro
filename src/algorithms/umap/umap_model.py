"""
Scikit-Learn Compatible UMAP Estimator implementation in Pure NumPy & SciPy with Multiprocessing support.
Provides fit, fit_transform, and out-of-sample transform.
"""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin

from .metrics import pairwise_distances, find_knn, resolve_n_jobs
from .fuzzy_simplicial_set import smooth_knn_dist, fuzzy_simplicial_set
from .curve_fitting import find_ab_params
from .layout import spectral_layout, optimize_layout

class UMAP(BaseEstimator, TransformerMixin):
    """
    Uniform Manifold Approximation and Projection (UMAP) for Dimension Reduction.
    
    Parameters
    ----------
    n_neighbors : int, default=15
        Number of nearest neighbors to construct fuzzy simplicial set graph.
    n_components : int, default=2
        Dimension of target embedded space.
    metric : str, default='euclidean'
        Distance metric ('euclidean', 'manhattan', 'cosine', 'tanimoto', 'jaccard', 'hamming', 'dice').
    min_dist : float, default=0.1
        Minimum distance between points in low-dimensional space.
    spread : float, default=1.0
        Effective spatial scale of points in low-dimensional space.
    n_epochs : int, default=100
        Number of optimization epochs. Reduced from 200 for faster convergence.
    learning_rate : float, default=1.0
        Initial learning rate for SGD.
    negative_sample_rate : int, default=5
        Number of negative samples per positive edge during optimization.
    init : str, default='spectral'
        Initialization method ('spectral', 'pca', or 'random').
    n_jobs : int, optional, default=None
        Number of parallel jobs/threads. If None or <=0, defaults to 50% of available CPU logical cores.
    random_state : int, default=42
        Random seed for reproducibility.
    """

    def __init__(
        self,
        n_neighbors: int = 15,
        n_components: int = 2,
        metric: str = 'euclidean',
        min_dist: float = 0.1,
        spread: float = 1.0,
        n_epochs: int = 100,
        learning_rate: float = 1.0,
        negative_sample_rate: int = 5,
        init: str = 'spectral',
        n_jobs: int = None,
        random_state: int = 42
    ):
        self.n_neighbors = n_neighbors
        self.n_components = n_components
        self.metric = metric
        self.min_dist = min_dist
        self.spread = spread
        self.n_epochs = n_epochs
        self.learning_rate = learning_rate
        self.negative_sample_rate = negative_sample_rate
        self.init = init
        self.n_jobs = n_jobs
        self.random_state = random_state

        # Fitted attributes
        self.embedding_ = None
        self.X_fit_ = None
        self.a_ = None
        self.b_ = None
        self.graph_ = None
        self.n_jobs_resolved_ = None

    def fit(self, X: np.ndarray, y=None):
        """Fits UMAP model to X."""
        self.fit_transform(X, y)
        return self

    def fit_transform(self, X: np.ndarray, y=None) -> np.ndarray:
        """Fits UMAP model to X and returns low-dimensional embeddings."""
        X = np.asarray(X, dtype=np.float64)
        self.X_fit_ = X
        self.n_jobs_resolved_ = resolve_n_jobs(self.n_jobs)

        n_samples = X.shape[0]
        if n_samples <= 1:
            self.embedding_ = np.zeros((n_samples, self.n_components))
            return self.embedding_

        effective_k = min(self.n_neighbors, n_samples - 1)

        # 1. Compute Pairwise Distances & KNN (Parallel)
        D = pairwise_distances(X, metric=self.metric, n_jobs=self.n_jobs_resolved_)
        knn_indices, knn_dists = find_knn(D, n_neighbors=effective_k, n_jobs=self.n_jobs_resolved_)

        # 2. Build Fuzzy Simplicial Set Graph (Parallel)
        self.graph_ = fuzzy_simplicial_set(
            D=D,
            knn_indices=knn_indices,
            knn_dists=knn_dists,
            n_neighbors=effective_k,
            n_jobs=self.n_jobs_resolved_
        )

        # 3. Fit (a, b) parameters
        self.a_, self.b_ = find_ab_params(spread=self.spread, min_dist=self.min_dist)

        # 4. Initialization
        if self.init == 'spectral':
            init_Y = spectral_layout(self.graph_, n_components=self.n_components, random_state=self.random_state)
        elif self.init == 'pca':
            from sklearn.decomposition import PCA
            init_Y = PCA(n_components=self.n_components, random_state=self.random_state).fit_transform(X) * 0.0001
        else:
            rng = np.random.default_rng(self.random_state)
            init_Y = 0.0001 * rng.standard_normal((n_samples, self.n_components))

        # 5. Optimize SGD Layout
        self.embedding_ = optimize_layout(
            graph=self.graph_,
            n_components=self.n_components,
            n_epochs=self.n_epochs,
            learning_rate=self.learning_rate,
            a=self.a_,
            b=self.b_,
            negative_sample_rate=self.negative_sample_rate,
            init_embedding=init_Y,
            random_state=self.random_state
        )

        return self.embedding_

    def transform(self, X_new: np.ndarray) -> np.ndarray:
        """Projects new unseen samples X_new into existing fitted UMAP embedding space."""
        if self.embedding_ is None or self.X_fit_ is None:
            raise ValueError("UMAP model is not fitted yet. Call fit_transform first.")

        X_new = np.asarray(X_new, dtype=np.float64)
        n_new = X_new.shape[0]

        # Distances between new points and fitted training points (Parallel)
        D_new = pairwise_distances(X_new, self.X_fit_, metric=self.metric, n_jobs=self.n_jobs_resolved_)
        k = min(self.n_neighbors, self.X_fit_.shape[0])

        Y_new = np.zeros((n_new, self.n_components), dtype=np.float64)

        for i in range(n_new):
            row = D_new[i]
            nn_idx = np.argsort(row)[:k]
            nn_dists = row[nn_idx]
            
            # Initial position is weighted average of nearest neighbors' embeddings
            weights = np.exp(-nn_dists / np.maximum(np.mean(nn_dists), 1e-12))
            weights_sum = np.sum(weights)
            if weights_sum > 0:
                weights /= weights_sum
                Y_new[i] = np.dot(weights, self.embedding_[nn_idx])
            else:
                Y_new[i] = np.mean(self.embedding_[nn_idx], axis=0)

            # Optimization for out-of-sample embedding
            alpha = self.learning_rate
            for epoch in range(30):
                for j_idx, neighbor_idx in enumerate(nn_idx):
                    w = weights[j_idx]
                    diff = Y_new[i] - self.embedding_[neighbor_idx]
                    dist2 = np.dot(diff, diff) + 1e-12

                    attr_coeff = (-2.0 * self.a_ * self.b_ * (dist2 ** (self.b_ - 1.0))) / (1.0 + self.a_ * (dist2 ** self.b_))
                    grad = np.clip(w * attr_coeff * diff, -4.0, 4.0)

                    Y_new[i] += alpha * grad

        return Y_new
