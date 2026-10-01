import numpy as np
from .base import BaseSplitter


class RandomSplitter(BaseSplitter):
    """Splits the data randomly to match the target training set ratio."""

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        n = len(df)
        indices = np.arange(n)
        np.random.shuffle(indices)
        split_idx = int(n * self.target_ratio)
        train_idx = indices[:split_idx]
        test_idx = indices[split_idx:]
        return list(train_idx), list(test_idx)
