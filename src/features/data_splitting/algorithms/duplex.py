import numpy as np
from .base import BaseSplitter
from ..utils.descriptors import preprocess_descriptors_numpy
from scipy.spatial.distance import cdist


class DuplexSplitter(BaseSplitter):
    """Alternates assigning furthest-apart points to the training and
    test sets, generalized to support arbitrary user-defined partition
    ratios."""

    def __init__(self, target_ratio: float, n_jobs: int = 1,
                 metric: str = 'euclidean', **kwargs):
        super().__init__(target_ratio, n_jobs, **kwargs)
        self.metric = metric

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        if descriptor_df is None or descriptor_df.empty:
            raise ValueError("Duplex splitting requires calculated molecular descriptors.")

        if self.metric == 'jaccard':
            X = preprocess_descriptors_numpy(descriptor_df, standardize=False).astype(bool)
        else:
            X = preprocess_descriptors_numpy(descriptor_df, standardize=True)

        n = X.shape[0]
        k = int(n * self.target_ratio)

        if k <= 0 or k >= n:
            return list(range(k)), list(range(k, n))

        dist_matrix = cdist(X, X, metric=self.metric)

        # Initialize training set with the two furthest points
        i1, j1 = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)
        train_idx = [i1, j1]
        remaining = set(range(n)) - {i1, j1}

        # Initialize test set with the two furthest remaining points
        if len(remaining) >= 2:
            rem_list = list(remaining)
            sub_dists = dist_matrix[rem_list][:, rem_list]
            r_i, r_j = np.unravel_index(np.argmax(sub_dists), sub_dists.shape)
            i2 = rem_list[r_i]
            j2 = rem_list[r_j]
            test_idx = [i2, j2]
            remaining -= {i2, j2}
        else:
            test_idx = list(remaining)
            remaining = set()

        # Distance tracks
        min_dist_train = (np.min(dist_matrix[:, train_idx], axis=1)
                          if len(train_idx) > 0 else np.zeros(n))
        min_dist_test = (np.min(dist_matrix[:, test_idx], axis=1)
                         if len(test_idx) > 0 else np.zeros(n))

        # Sequentially assign remaining points to maintain the target ratio
        while remaining:
            current_ratio = len(train_idx) / (len(train_idx) + len(test_idx))
            rem_list = list(remaining)

            if current_ratio < self.target_ratio and len(train_idx) < k:
                # Add to train: select point furthest from current training set
                best_idx = rem_list[np.argmax(min_dist_train[rem_list])]
                train_idx.append(best_idx)
                remaining.remove(best_idx)
                min_dist_train = np.minimum(min_dist_train, dist_matrix[:, best_idx])
            else:
                # Add to test: select point furthest from current test set
                best_idx = rem_list[np.argmax(min_dist_test[rem_list])]
                test_idx.append(best_idx)
                remaining.remove(best_idx)
                min_dist_test = np.minimum(min_dist_test, dist_matrix[:, best_idx])

        return train_idx, test_idx
