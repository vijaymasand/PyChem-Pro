import pandas as pd


class BaseSplitter:
    """Base class for all QSAR data splitting algorithms."""

    def __init__(self, target_ratio: float, n_jobs: int = 1, **kwargs):
        self.target_ratio = target_ratio  # fraction of train set, e.g., 0.8
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def split(self, df: pd.DataFrame, smiles_col: str,
              target_col: str = None, descriptor_df: pd.DataFrame = None):
        """
        Perform splitting.
        Returns a tuple: (train_indices, test_indices)
        """
        raise NotImplementedError
