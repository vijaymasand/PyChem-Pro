import numpy as np
from .base import BaseSplitter
from ..utils.descriptors import preprocess_descriptors_numpy
from scipy.spatial.distance import cdist


class CADEXSplitter(BaseSplitter):
    """Uniformly covers the descriptor space by picking the two furthest
    points and iteratively adding points that maximize the minimum
    distance to already selected points (Kennard-Stone algorithm)."""

    def __init__(self, target_ratio: float, n_jobs: int = 1,
                 metric: str = 'euclidean', **kwargs):
        super().__init__(target_ratio, n_jobs, **kwargs)
        self.metric = metric

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        if descriptor_df is None or descriptor_df.empty:
            raise ValueError("CADEX (Kennard-Stone) requires calculated molecular descriptors.")

        if self.metric == 'jaccard':
            X = preprocess_descriptors_numpy(descriptor_df, standardize=False).astype(bool)
        else:
            X = preprocess_descriptors_numpy(descriptor_df, standardize=True)

        n = X.shape[0]
        k = int(n * self.target_ratio)

        if k <= 0 or k >= n:
            return list(range(k)), list(range(k, n))

        # Compute pairwise distances
        dist_matrix = cdist(X, X, metric=self.metric)

        # Select the two points with the maximum mutual distance
        i, j = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)

        selected = [i, j]
        remaining = set(range(n)) - {i, j}

        # Pre-calculate minimum distances to selected set
        min_distances = np.min(dist_matrix[:, selected], axis=1)

        # Iteratively select the furthest remaining point
        while len(selected) < k:
            rem_list = list(remaining)
            best_idx = rem_list[np.argmax(min_distances[rem_list])]

            selected.append(best_idx)
            remaining.remove(best_idx)

            # Update min_distances incrementally
            min_distances = np.minimum(min_distances, dist_matrix[:, best_idx])

        train_idx = selected
        test_idx = list(remaining)
        return train_idx, test_idx
