from dataclasses import dataclass
from typing import List, Dict, Any, Tuple, Optional
import pandas as pd
import numpy as np

from .algorithms import get_splitter, ALGOS_NEED_FEATURES, SPLITTERS
from .utils.descriptors import DescriptorMode, get_feature_matrix
from .utils.visualization import compute_split_stats, generate_pca_plot


@dataclass
class SplitResult:
    """Result of a data splitting operation."""
    train_indices: List[int]
    test_indices: List[int]
    annotated_df: pd.DataFrame
    stats: Dict[str, Any]
    algorithm: str
    target_ratio: float


class DataSplitEngine:
    """
    Headless engine for QSAR dataset splitting.
    This class provides a clean programmatic API that can be used from Python
    scripts or Jupyter notebooks, independent of any GUI.
    """

    def __init__(self):
        self.df: Optional[pd.DataFrame] = None
        self.original_columns: List[str] = []
        self._desc_df: Optional[pd.DataFrame] = None

    def load_csv(self, file_path: str) -> bool:
        """Load a CSV dataset into the engine."""
        try:
            self.df = pd.read_csv(file_path)
        except Exception:
            try:
                self.df = pd.read_csv(file_path, encoding='latin1')
            except Exception as e:
                raise ValueError(f"Failed to load CSV: {e}")
        
        self.original_columns = list(self.df.columns)
        return True

    def load_dataframe(self, df: pd.DataFrame) -> bool:
        """Load an existing pandas DataFrame into the engine."""
        self.df = df.copy()
        self.original_columns = list(self.df.columns)
        return True

    def guess_columns(self) -> Dict[str, Optional[str]]:
        """Attempt to auto-detect SMILES, ID, and target columns."""
        if self.df is None:
            return {"smiles": None, "id": None, "target": None}
            
        cols = [c.lower() for c in self.original_columns]
        result = {"smiles": None, "id": None, "target": None}
        
        # Guess SMILES
        for i, c in enumerate(cols):
            if "smiles" in c:
                result["smiles"] = self.original_columns[i]
                break
                
        # Guess ID
        for i, c in enumerate(cols):
            if any(k in c for k in ("id", "name", "title", "compound")):
                result["id"] = self.original_columns[i]
                break
                
        # Guess Target
        for i, c in enumerate(cols):
            if any(k in c for k in ("activity", "target", "value", "pic50", "ic50", "pki", "pec50", "logd", "logp")):
                result["target"] = self.original_columns[i]
                break
                
        return result

    def get_numeric_columns(self) -> List[str]:
        """Return a list of all numeric columns in the dataset."""
        if self.df is None:
            return []
        numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
        return list(self.df.select_dtypes(include=numerics).columns)

    @staticmethod
    def available_algorithms() -> List[str]:
        """List all available splitting algorithms."""
        return list(SPLITTERS.keys())

    def split(self,
              algorithm: str,
              target_ratio: float,
              smiles_col: str,
              id_col: Optional[str] = None,
              target_col: Optional[str] = None,
              desc_mode: DescriptorMode = DescriptorMode.FINGERPRINTS,
              csv_desc_cols: Optional[List[str]] = None,
              fp_type: Optional[str] = None,
              metric: str = "euclidean",
              n_jobs: int = 1) -> SplitResult:
        """
        Perform dataset splitting.

        Parameters
        ----------
        algorithm : str
            Name of the algorithm (e.g. "Random", "Sphere Exclusion").
        target_ratio : float
            Fraction of data to assign to the training set (e.g. 0.8).
        smiles_col : str
            Column name containing SMILES strings.
        id_col : str, optional
            Column name containing molecule IDs.
        target_col : str, optional
            Column name containing target/activity values.
        desc_mode : DescriptorMode
            How to obtain features for distance-based algorithms.
        csv_desc_cols : list of str, optional
            Columns to use if desc_mode == USE_CSV.
        fp_type : str, optional
            Fingerprint type to use if desc_mode == FINGERPRINTS.
        metric : str
            Distance metric to use ('euclidean', 'jaccard', etc).
        n_jobs : int
            Number of parallel workers.

        Returns
        -------
        SplitResult
            Dataclass containing indices, stats, and the annotated DataFrame.
        """
        if self.df is None:
            raise ValueError("No data loaded. Call load_csv() or load_dataframe() first.")
            
        if algorithm not in SPLITTERS:
            raise ValueError(f"Unknown algorithm: {algorithm}. Available: {self.available_algorithms()}")
            
        if smiles_col not in self.df.columns:
            raise ValueError(f"SMILES column '{smiles_col}' not found in dataset.")

        # Data cleansing
        df_clean = self.df.dropna(subset=[smiles_col]).copy()
        df_clean[smiles_col] = df_clean[smiles_col].astype(str).str.strip()
        df_clean = df_clean[df_clean[smiles_col] != ""]
        
        initial_len = len(df_clean)
        df_clean = df_clean.drop_duplicates(subset=[smiles_col]).reset_index(drop=True)
        duplicates_dropped = initial_len - len(df_clean)
        
        if len(df_clean) == 0:
            raise ValueError("Dataset contains no valid SMILES strings after cleansing.")

        smiles_list = df_clean[smiles_col].tolist()
        names_list = df_clean[id_col].tolist() if id_col and id_col in df_clean.columns else [f"Mol_{i}" for i in range(len(df_clean))]

        needs_features = algorithm in ALGOS_NEED_FEATURES
        self._desc_df = pd.DataFrame()
        valid_indices = list(range(len(df_clean)))

        if needs_features:
            self._desc_df, valid_indices = get_feature_matrix(
                df_clean, smiles_list, names_list,
                mode=desc_mode,
                csv_descriptor_cols=csv_desc_cols,
                fp_type=fp_type,
                n_jobs=n_jobs
            )

            if self._desc_df is None or self._desc_df.empty:
                raise RuntimeError("Feature matrix is empty. Check SMILES / descriptor columns.")

            df_clean = df_clean.iloc[valid_indices].reset_index(drop=True)
            self._desc_df = self._desc_df.reset_index(drop=True)

        # Initialize and run splitter
        splitter = get_splitter(algorithm, target_ratio=target_ratio, n_jobs=n_jobs, metric=metric)
        train_idx, test_idx = splitter.split(
            df_clean, smiles_col=smiles_col,
            target_col=target_col, descriptor_df=self._desc_df
        )

        # Annotate dataframe
        status = ["Test"] * len(df_clean)
        for i in train_idx:
            status[i] = "Train"
        df_clean["Split_Status"] = status

        # Compute statistics
        stats = compute_split_stats(df_clean, train_idx, test_idx, target_col=target_col)
        stats["duplicates_dropped"] = duplicates_dropped

        return SplitResult(
            train_indices=train_idx,
            test_indices=test_idx,
            annotated_df=df_clean,
            stats=stats,
            algorithm=algorithm,
            target_ratio=target_ratio
        )

    def plot_pca(self, train_indices: List[int], test_indices: List[int], save_path: str) -> bool:
        """
        Generate a PCA scatter plot of the split.
        Must be called after `split()` if the algorithm required features.
        If the algorithm didn't require features, `split()` wouldn't have computed them.
        """
        if self._desc_df is None or self._desc_df.empty:
            raise RuntimeError("No feature matrix available to perform PCA. Was a distance-based splitting algorithm used?")
            
        return generate_pca_plot(self._desc_df, train_indices, test_indices, save_path)


def split_dataset(df_or_path, algorithm: str, target_ratio: float, smiles_col: str, **kwargs) -> SplitResult:
    """
    Convenience function for one-shot dataset splitting.
    
    Example:
        result = split_dataset("data.csv", "Sphere Exclusion", 0.8, "SMILES", target_col="pIC50")
        result.annotated_df.to_csv("split_data.csv")
    """
    engine = DataSplitEngine()
    if isinstance(df_or_path, str):
        engine.load_csv(df_or_path)
    elif isinstance(df_or_path, pd.DataFrame):
        engine.load_dataframe(df_or_path)
    else:
        raise TypeError("df_or_path must be a file path string or a pandas DataFrame.")
        
    return engine.split(algorithm, target_ratio, smiles_col, **kwargs)
