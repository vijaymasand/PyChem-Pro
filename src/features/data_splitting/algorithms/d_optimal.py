import numpy as np
from .base import BaseSplitter
from ..utils.descriptors import preprocess_descriptors_numpy


class DOptimalSplitter(BaseSplitter):
    """Greedily selects training samples to maximize the determinant of
    the information matrix X^T X, ensuring high statistical efficiency."""

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        if descriptor_df is None or descriptor_df.empty:
            raise ValueError("D-optimal splitting requires calculated molecular descriptors.")

        X = preprocess_descriptors_numpy(descriptor_df)
        n = X.shape[0]
        k = int(n * self.target_ratio)

        if k <= 0 or k >= n:
            return list(range(k)), list(range(k, n))

        d = X.shape[1]
        lam = 1e-4  # Ridge regularization parameter to ensure non-singularity

        # Start with the point furthest from the origin (max norm)
        norms = np.linalg.norm(X, axis=1)
        first_idx = np.argmax(norms)

        selected = [first_idx]
        remaining = set(range(n)) - {first_idx}

        # Sequentially add points that maximize prediction variance (leverage)
        while len(selected) < k:
            X_S = X[selected]
            # Compute inverse of information matrix
            info_matrix = np.dot(X_S.T, X_S) + lam * np.eye(d)
            try:
                inv_info = np.linalg.inv(info_matrix)
            except np.linalg.LinAlgError:
                inv_info = np.linalg.pinv(info_matrix)

            rem_list = list(remaining)
            X_rem = X[rem_list]

            # Vectorized score: s_i = x_i^T * inv_info * x_i
            scores = np.sum((X_rem @ inv_info) * X_rem, axis=1)

            best_local_idx = np.argmax(scores)
            best_global_idx = rem_list[best_local_idx]

            selected.append(best_global_idx)
            remaining.remove(best_global_idx)

        train_idx = selected
        test_idx = list(remaining)
        return train_idx, test_idx
