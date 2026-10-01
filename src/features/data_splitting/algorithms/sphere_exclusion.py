import numpy as np
from .base import BaseSplitter
from ..utils.descriptors import preprocess_descriptors_numpy
from scipy.spatial.distance import cdist


def run_sphere_exclusion(X, R, metric='euclidean'):
    """Runs a single pass of the sphere exclusion algorithm with radius R."""
    n = X.shape[0]
    dist_matrix = cdist(X, X, metric=metric)

    pool = set(range(n))
    train_idx = []
    test_idx = []

    while pool:
        pool_list = list(pool)
        # Select the seed compound (closest to the centroid of the remaining pool)
        if len(pool_list) > 1:
            sub_dists = dist_matrix[pool_list][:, pool_list]
            seed_idx = pool_list[np.argmin(np.mean(sub_dists, axis=1))]
        else:
            seed_idx = pool_list[0]

        train_idx.append(seed_idx)
        pool.remove(seed_idx)

        # Exclude neighbors within distance R and assign them to the test set
        neighbors = [p for p in pool if dist_matrix[seed_idx, p] <= R]
        for neighbor in neighbors:
            test_idx.append(neighbor)
            pool.remove(neighbor)

    return train_idx, test_idx


class SphereExclusionSplitter(BaseSplitter):
    """Performs a binary search over the exclusion sphere radius R to
    find the configuration that matches the target ratio as closely
    as possible."""

    def __init__(self, target_ratio: float, n_jobs: int = 1,
                 metric: str = 'euclidean', **kwargs):
        super().__init__(target_ratio, n_jobs, **kwargs)
        self.metric = metric

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        if descriptor_df is None or descriptor_df.empty:
            raise ValueError("Sphere Exclusion splitting requires calculated molecular descriptors.")

        if self.metric == 'jaccard':
            X = preprocess_descriptors_numpy(descriptor_df, standardize=False).astype(bool)
        else:
            X = preprocess_descriptors_numpy(descriptor_df, standardize=True)

        n = X.shape[0]
        target_train_size = int(n * self.target_ratio)

        # Compute distances to find max distance for search bounds
        dist_matrix = cdist(X, X, metric=self.metric)
        max_dist = np.max(dist_matrix)

        # Binary search for the optimal exclusion radius R
        low = 0.0
        high = max_dist
        best_R = 0.0
        best_train_idx = list(range(n))
        best_test_idx = []
        best_diff = n

        # Perform 15 iterations of binary search
        for _ in range(15):
            mid = (low + high) / 2.0
            train_idx, test_idx = run_sphere_exclusion(X, mid, metric=self.metric)

            diff = len(train_idx) - target_train_size
            if abs(diff) < best_diff:
                best_diff = abs(diff)
                best_R = mid
                best_train_idx = train_idx
                best_test_idx = test_idx

            # If training set is too large, we need a larger radius to exclude more points
            if len(train_idx) > target_train_size:
                low = mid
            else:
                high = mid

        return best_train_idx, best_test_idx
